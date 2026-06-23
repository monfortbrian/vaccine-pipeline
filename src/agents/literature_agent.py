"""
LITERATURE AGENT
Searches published evidence for each predicted epitope candidate.

Tools:
  PubMed E-utilities API  : fetch abstracts by protein/epitope query
  Qdrant (in-memory)      : vector store for semantic similarity search
  sentence-transformers   : embed abstracts and epitope sequences
  Claude API (optional)   : synthesize evidence into structured summary
"""

import os
import re
import time
import logging
import hashlib
import requests
from typing import List, Dict, Any, Optional

from src.models.candidate import CandidateProtein, ConfidenceTier
from src.utils.logger import get_logger

logger = get_logger("tope_deep.agents.Agent 9")

PUBMED_API_KEY = os.getenv("NCBI_API_KEY", "")
PUBMED_BASE    = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
MAX_ABSTRACTS  = 20
EMBED_MODEL    = "all-MiniLM-L6-v2"
QDRANT_URL     = os.getenv("QDRANT_URL", "")


# ── Query builder ─────────────────────────────────────────────────────────────

_VACCINE_TERMS = (
    '"T cell"[Title/Abstract] OR '
    '"epitope"[Title/Abstract] OR '
    '"vaccine"[Title/Abstract] OR '
    '"immunogen"[Title/Abstract] OR '
    '"antigen"[Title/Abstract] OR '
    '"MHC"[Title/Abstract] OR '
    '"antibody"[Title/Abstract]'
)

_UNIPROT_ID_RE = re.compile(r'^[A-Z][0-9][A-Z0-9]{3}[0-9]$|^[OPQ][0-9][A-Z0-9]{3}[0-9]$')


def _clean_protein_name(name: str) -> str:
    """
    Strip UniProt parenthetical accession from protein name before PubMed query.
    "Mycolipanoate synthase (A0A089QRB9)" to "Mycolipanoate synthase"
    PubMed does not index UniProt accessions in freetext including them
    causes zero results.
    """
    return re.sub(r'\s*\([A-Z0-9]{6,12}\)\s*$', '', name).strip()


def _is_real_uniprot_id(protein_id: str) -> bool:
    """True if protein_id looks like a UniProt accession, not 'user_input'."""
    if not protein_id or protein_id == "user_input":
        return False
    return bool(re.match(r'^[A-Z0-9]{6,10}$', protein_id))


def _build_pubmed_query(candidate) -> str:
    """
    Build a PubMed query adaptive to any input type.
    Never hardcodes organism or pathogen name.
    Always strips UniProt parenthetical from protein name.
    """
    protein_name = candidate.protein_name or ""
    protein_id   = candidate.protein_id   or ""
    source       = getattr(candidate, 'source', 'uniprot')
    decisions    = getattr(candidate, 'decisions', [])

    clean_name  = _clean_protein_name(protein_name)
    has_uniprot = _is_real_uniprot_id(protein_id)

    # Case 1: user-submitted sequence with generic name
    if source == "user_input" or protein_id == "user_input":
        if not clean_name or clean_name.lower() in ("custom protein", "unknown", "protein"):
            ctl_epitopes = [
                ep.sequence for ep in (candidate.ctl_epitopes or [])
                if hasattr(ep, 'sequence')
            ][:2]
            if ctl_epitopes:
                ep_terms = " OR ".join(f'"{seq}"[Title/Abstract]' for seq in ctl_epitopes)
                return f'({ep_terms}) AND ({_VACCINE_TERMS})'
            return f'("vaccine epitope"[Title/Abstract]) AND ({_VACCINE_TERMS})'

    # Case 2: real protein name, with or without UniProt ID
    name_term = f'"{clean_name}"[Title/Abstract]'

    if has_uniprot:
        id_term = f'"{protein_id}"[Title/Abstract]'
        subject = f'({name_term} OR {id_term})'
    else:
        subject = f'({name_term})'

    return f'{subject} AND ({_VACCINE_TERMS})'


def _build_fallback_query(candidate) -> str:
    """
    Fallback when primary returns 0 results.
    Drops vaccine terms, searches name OR accession alone.
    """
    protein_name = candidate.protein_name or ""
    protein_id   = candidate.protein_id   or ""

    clean_name  = _clean_protein_name(protein_name)
    has_uniprot = _is_real_uniprot_id(protein_id)

    if has_uniprot and clean_name:
        return f'("{clean_name}"[Title/Abstract] OR "{protein_id}"[Title/Abstract])'
    elif has_uniprot:
        return f'"{protein_id}"[Title/Abstract]'
    elif clean_name and clean_name.lower() not in ("custom protein", "unknown"):
        return f'"{clean_name}"[Title/Abstract]'
    else:
        return ""


# ── PubMed fetch ──────────────────────────────────────────────────────────────

def _fetch_pubmed_abstracts(query: str, max_results: int = MAX_ABSTRACTS) -> List[Dict]:
    try:
        search_resp = requests.get(
            f"{PUBMED_BASE}/esearch.fcgi",
            params={
                "db":      "pubmed",
                "term":    query,
                "retmax":  max_results,
                "retmode": "json",
                "api_key": PUBMED_API_KEY,
            },
            timeout=15,
        )
        search_resp.raise_for_status()
        pmids = search_resp.json().get("esearchresult", {}).get("idlist", [])

        if not pmids:
            logger.info(f"Agent 9: PubMed returned 0 PMIDs for query: {query[:80]}")
            return []

        time.sleep(0.1 if PUBMED_API_KEY else 0.34)

        fetch_resp = requests.get(
            f"{PUBMED_BASE}/efetch.fcgi",
            params={
                "db":      "pubmed",
                "id":      ",".join(pmids),
                "rettype": "abstract",
                "retmode": "xml",
                "api_key": PUBMED_API_KEY,
            },
            timeout=20,
        )
        fetch_resp.raise_for_status()

        import xml.etree.ElementTree as ET
        root = ET.fromstring(fetch_resp.content)
        articles = []
        for article in root.findall(".//PubmedArticle"):
            pmid_el      = article.find(".//PMID")
            title_el     = article.find(".//ArticleTitle")
            abstract_els = article.findall(".//AbstractText")
            year_el      = article.find(".//PubDate/Year")
            journal_el   = article.find(".//Journal/Title")

            pmid     = pmid_el.text     if pmid_el     is not None else ""
            title    = title_el.text    if title_el    is not None else ""
            abstract = " ".join((el.text or "") for el in abstract_els).strip()
            year     = year_el.text     if year_el     is not None else ""
            journal  = journal_el.text  if journal_el  is not None else ""

            if abstract:
                articles.append({
                    "pmid": pmid, "title": title,
                    "abstract": abstract, "journal": journal, "year": year,
                })

        logger.info(f"Agent 9: PubMed returned {len(articles)} abstracts for: {query[:60]}")
        return articles

    except Exception as e:
        logger.warning(f"Agent 9: PubMed fetch failed: {e}")
        return []


# ── Vector store ──────────────────────────────────────────────────────────────

def _get_qdrant_client():
    try:
        from qdrant_client import QdrantClient
        if QDRANT_URL:
            logger.info(f"Agent 9: Qdrant {QDRANT_URL}")
            return QdrantClient(url=QDRANT_URL)
        logger.info("Agent 9: Qdrant in-memory")
        return QdrantClient(":memory:")
    except ImportError:
        logger.warning("Agent 9: qdrant-client not installed, falling back to ChromaDB")
        return None


def _get_chroma_client():
    try:
        import chromadb
        return chromadb.Client()
    except ImportError:
        logger.error("Agent 9: Neither qdrant-client nor chromadb installed")
        return None


def _get_embedder():
    try:
        from sentence_transformers import SentenceTransformer
        return SentenceTransformer(EMBED_MODEL)
    except ImportError:
        logger.warning("Agent 9: sentence-transformers not installed")
        return None


def _build_collection_name(run_id: str, protein_id: str) -> str:
    h = hashlib.md5(f"{run_id}_{protein_id}".encode()).hexdigest()[:8]
    return f"td_{h}"


def _index_abstracts_qdrant(client, collection_name: str, abstracts: List[Dict], embedder) -> bool:
    try:
        from qdrant_client.models import Distance, VectorParams, PointStruct
        texts   = [f"{a['title']} {a['abstract']}" for a in abstracts]
        vectors = embedder.encode(texts, show_progress_bar=False).tolist()
        client.recreate_collection(
            collection_name=collection_name,
            vectors_config=VectorParams(size=len(vectors[0]), distance=Distance.COSINE),
        )
        client.upsert(
            collection_name=collection_name,
            points=[
                PointStruct(id=i, vector=vectors[i], payload=abstracts[i])
                for i in range(len(abstracts))
            ],
        )
        return True
    except Exception as e:
        logger.warning(f"Agent 9: Qdrant index failed: {e}")
        return False


def _search_evidence(client, collection_name: str, query_text: str, embedder, top_k: int = 3) -> List[Dict]:
    try:
        query_vector = embedder.encode([query_text], show_progress_bar=False)[0].tolist()
        results = client.search(
            collection_name=collection_name,
            query_vector=query_vector,
            limit=top_k,
            score_threshold=0.3,
        )
        return [r.payload for r in results]
    except Exception as e:
        logger.warning(f"Agent 9: Qdrant search failed: {e}")
        return []


def _detect_failure_signals(abstracts: List[Dict]) -> List[str]:
    failure_keywords = {
        "immune evasion":     "immune evasion reported in literature",
        "highly variable":    "high sequence variability escape risk",
        "polymorphic":        "polymorphic region population coverage may vary",
        "cross-reactive":     "cross-reactivity with host proteins reported",
        "poor immunogen":     "poor immunogenicity in human trials",
        "no t cell response": "absence of T-cell response in human subjects",
        "failed phase":       "clinical trial failure reported",
        "toxic":              "toxicity concerns noted",
        "allergen":           "allergenicity signal in literature",
    }
    found = []
    for abstract in abstracts:
        text = (abstract.get("abstract", "") + " " + abstract.get("title", "")).lower()
        for keyword, signal in failure_keywords.items():
            if keyword in text and signal not in found:
                found.append(signal)
    return found


def _synthesize_with_claude(protein_name: str, abstracts: List[Dict], failure_signals: List[str]) -> Optional[str]:
    api_key = os.getenv("ANTHROPIC_API_KEY", "")
    if not api_key or not abstracts:
        return None
    try:
        abstract_text = "\n\n".join([
            f"[PMID {a['pmid']}] {a['title']}\n{a['abstract'][:400]}"
            for a in abstracts[:5]
        ])
        failure_text = (
            f"Known failure signals: {', '.join(failure_signals)}"
            if failure_signals else "No failure signals detected."
        )
        prompt = f"""You are a computational immunologist reviewing published literature
for early vaccine discovery. Summarize the evidence below for {protein_name}
in 2-3 sentences. State: (1) what prior evidence exists for this protein as a vaccine target,
(2) any known concerns or failure signals, (3) your assessment of evidence quality.
Be precise and scientific. No hedging language.

{failure_text}

Literature:
{abstract_text}"""

        resp = requests.post(
            "https://api.anthropic.com/v1/messages",
            headers={
                "x-api-key":         api_key,
                "anthropic-version": "2023-06-01",
                "content-type":      "application/json",
            },
            json={
                "model":    "claude-sonnet-4-6",
                "max_tokens": 300,
                "messages": [{"role": "user", "content": prompt}],
            },
            timeout=20,
        )
        resp.raise_for_status()
        return resp.json()["content"][0]["text"].strip()
    except Exception as e:
        logger.warning(f"Agent 9: Claude synthesis failed: {e}")
        return None


# ── Agent ─────────────────────────────────────────────────────────────────────

class LiteratureAgent:
    def __init__(self):
        self.stage_name  = "literature_search"
        self._qdrant     = None
        self._embedder   = None
        self._use_qdrant = True

    def _init_clients(self):
        if self._qdrant is None:
            self._qdrant = _get_qdrant_client()
            if self._qdrant is None:
                self._qdrant     = _get_chroma_client()
                self._use_qdrant = False
        if self._embedder is None:
            self._embedder = _get_embedder()

    def run(self, candidates: List[CandidateProtein], run_id: str = "unknown") -> List[CandidateProtein]:
        logger.info("Agent 9: Starting literature search")
        self._init_clients()

        active = [c for c in candidates if c.status.value == "active"]
        logger.info(f"   {len(active)} candidates")

        for i, candidate in enumerate(active):
            logger.info(f"   [{i+1}/{len(active)}] {candidate.protein_name}")
            start = time.time()

            # FIX: define clean_name BEFORE it is used anywhere in this loop
            clean_name = _clean_protein_name(candidate.protein_name or "")

            # ── PRIMARY QUERY ────────────────────────────────────────────────
            primary_query = _build_pubmed_query(candidate)
            logger.info(f"      Query: {primary_query[:120]}")
            abstracts = _fetch_pubmed_abstracts(primary_query)

            # ── FALLBACK QUERY (if 0 results) ─────────────────────────────
            if not abstracts:
                fallback_query = _build_fallback_query(candidate)
                if fallback_query:
                    logger.info(f"      Fallback query: {fallback_query[:120]}")
                    abstracts = _fetch_pubmed_abstracts(fallback_query, max_results=10)
                else:
                    logger.info("      No fallback query available (user sequence with no name)")

            # ── EPITOPE-LEVEL QUERIES ─────────────────────────────────────
            high_conf_epitopes = [
                ep for ep in candidate.ctl_epitopes
                if ep.confidence_tier == ConfidenceTier.HIGH
            ][:3]

            epitope_abstracts = []
            for ep in high_conf_epitopes:
                ep_query = f'"{ep.sequence}"[Title/Abstract] AND "T cell"[Title/Abstract]'
                ep_abs   = _fetch_pubmed_abstracts(ep_query, max_results=5)
                epitope_abstracts.extend(ep_abs)

            # ── DEDUPLICATE ───────────────────────────────────────────────
            all_abstracts    = abstracts + epitope_abstracts
            seen_pmids       = set()
            unique_abstracts = []
            for a in all_abstracts:
                if a["pmid"] not in seen_pmids:
                    seen_pmids.add(a["pmid"])
                    unique_abstracts.append(a)

            logger.info(f"      {len(unique_abstracts)} unique abstracts")

            failure_signals = _detect_failure_signals(unique_abstracts)
            prior_validated = any(
                ep.sequence in (a.get("abstract", "") + a.get("title", ""))
                for ep in high_conf_epitopes
                for a in unique_abstracts
            )

            # ── INDEX IN QDRANT ───────────────────────────────────────────
            collection_name = _build_collection_name(run_id, candidate.protein_id)
            if self._qdrant and self._embedder and unique_abstracts and self._use_qdrant:
                indexed = _index_abstracts_qdrant(
                    self._qdrant, collection_name, unique_abstracts, self._embedder
                )
                if indexed:
                    logger.info(f"      Indexed {len(unique_abstracts)} abstracts in Qdrant")

            # ── CLAUDE SYNTHESIS ──────────────────────────────────────────
            # FIX: clean_name is now defined above, no longer used before assignment
            literature_summary = _synthesize_with_claude(
                clean_name,
                unique_abstracts,
                failure_signals,
            )

            elapsed = round(time.time() - start, 1)

            reasoning_parts = [
                f"PubMed search for '{clean_name}' returned {len(unique_abstracts)} abstracts.",
                f"Prior experimental validation: {prior_validated}.",
            ]
            if failure_signals:
                reasoning_parts.append(f"Failure signals: {'; '.join(failure_signals)}.")
            if unique_abstracts:
                pmid_list = ', '.join(a['pmid'] for a in unique_abstracts[:5])
                reasoning_parts.append(f"Evidence PMIDs: {pmid_list}.")
            if literature_summary:
                reasoning_parts.append(literature_summary)

            candidate.add_decision(
                stage=self.stage_name,
                decision="literature_searched",
                reasoning=" ".join(reasoning_parts),
                pubmed_hits=len(unique_abstracts),
                prior_validated=prior_validated,
                evidence_pmids=[a["pmid"] for a in unique_abstracts[:10]],
                failure_signals=failure_signals,
                literature_summary=literature_summary,
                search_time_s=elapsed,
                query=primary_query,
                result_count=len(unique_abstracts),
            )

            logger.info(
                f"      prior_validated={prior_validated} | "
                f"failure_signals={len(failure_signals)} | "
                f"pmids={len(unique_abstracts)} | {elapsed}s"
            )

        logger.info("Agent 9: Literature search complete")
        return candidates

    def get_status(self) -> Dict[str, Any]:
        return {
            "qdrant_available":   self._qdrant is not None,
            "embedder_available": self._embedder is not None,
            "qdrant_url":         QDRANT_URL or "in-memory",
            "pubmed_api_key":     bool(PUBMED_API_KEY),
        }


literature_agent = LiteratureAgent()