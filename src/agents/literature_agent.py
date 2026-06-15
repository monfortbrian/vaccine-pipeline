"""
LITERATURE AGENT
Searches published evidence for each predicted epitope candidate.

Architecture: retrieval-augmented search over PubMed abstracts.

Tools:
  PubMed E-utilities API  : fetch abstracts by protein/epitope query
  Qdrant (in-memory)      : vector store for semantic similarity search
  sentence-transformers   : embed abstracts and epitope sequences
  Claude API (optional)   : synthesize evidence signals into structured summary

Vector store strategy:
  - Run-scoped: each pipeline run gets its own Qdrant collection
  - In-memory by default (QdrantClient(":memory:"))
  - Persisted Qdrant service if QDRANT_URL env var is set
  - ChromaDB fallback if neither Qdrant option is available

For each candidate protein:
  1. Query PubMed for protein name + "epitope" + "T cell" OR "vaccine"
  2. Fetch abstracts (up to 20 per candidate)
  3. Embed abstracts with sentence-transformers (all-MiniLM-L6-v2)
  4. Index in Qdrant collection scoped to this run
  5. For each high-confidence CTL epitope, query for prior evidence
  6. Record: prior_validated (bool), failure_signals (list), evidence_pmids (list)
  7. Optional: Claude synthesizes evidence into one-paragraph summary per candidate

Output fields on each candidate decision record:
  pubmed_hits        : number of abstracts retrieved
  prior_validated    : whether any epitope appears in published assay data
  evidence_pmids     : list of supporting PMIDs
  failure_signals    : list of known failure modes found in literature
  literature_summary : Claude synthesis (if ANTHROPIC_API_KEY set)

References:
  PubMed E-utilities: Sayers et al. (2022) NCBI Insights
  sentence-transformers: Reimers & Gurevych (2019) EMNLP
  Qdrant: qdrant.tech (Apache 2.0)
"""

import os
import time
import logging
import hashlib
import requests
from typing import List, Dict, Any, Optional

from src.models.candidate import CandidateProtein, ConfidenceTier

from src.utils.logger import get_logger
logger = get_logger("tope_deep.agents.N10")  # use the correct agent name

PUBMED_API_KEY = os.getenv("NCBI_API_KEY", "")
PUBMED_BASE    = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
MAX_ABSTRACTS  = 20
EMBED_MODEL    = "all-MiniLM-L6-v2"   # 384-dim, fast, no GPU needed
QDRANT_URL     = os.getenv("QDRANT_URL", "")  # empty = in-memory


def _get_qdrant_client():
    """Return Qdrant client persisted if QDRANT_URL set, else in-memory."""
    try:
        from qdrant_client import QdrantClient
        if QDRANT_URL:
            logger.info(f"N9: Qdrant client -> {QDRANT_URL}")
            return QdrantClient(url=QDRANT_URL)
        logger.info("N9: Qdrant client -> in-memory")
        return QdrantClient(":memory:")
    except ImportError:
        logger.warning("N9: qdrant-client not installed, falling back to ChromaDB")
        return None


def _get_chroma_client():
    """ChromaDB fallback in-process, no separate service needed."""
    try:
        import chromadb
        logger.info("N9: ChromaDB client -> in-process fallback")
        return chromadb.Client()
    except ImportError:
        logger.error("N9: Neither qdrant-client nor chromadb installed")
        return None


def _get_embedder():
    """Load sentence-transformer model. Cached after first load."""
    try:
        from sentence_transformers import SentenceTransformer
        return SentenceTransformer(EMBED_MODEL)
    except ImportError:
        logger.warning("N9: sentence-transformers not installed, skipping semantic search")
        return None


def _fetch_pubmed_abstracts(query: str, max_results: int = MAX_ABSTRACTS) -> List[Dict]:
    """
    Search PubMed and fetch abstracts.
    Returns list of {pmid, title, abstract, journal, year}.
    """
    try:
        # Step 1: search for PMIDs
        search_resp = requests.get(
            f"{PUBMED_BASE}/esearch.fcgi",
            params={
                "db": "pubmed", "term": query,
                "retmax": max_results, "retmode": "json",
                "api_key": PUBMED_API_KEY,
            },
            timeout=15,
        )
        search_resp.raise_for_status()
        pmids = search_resp.json().get("esearchresult", {}).get("idlist", [])
        if not pmids:
            return []

        time.sleep(0.1 if PUBMED_API_KEY else 0.34)  # rate limit: 10/s with key, 3/s without

        # Step 2: fetch abstracts
        fetch_resp = requests.get(
            f"{PUBMED_BASE}/efetch.fcgi",
            params={
                "db": "pubmed", "id": ",".join(pmids),
                "rettype": "abstract", "retmode": "xml",
                "api_key": PUBMED_API_KEY,
            },
            timeout=20,
        )
        fetch_resp.raise_for_status()

        # Parse XML minimally extract title and abstract text
        import xml.etree.ElementTree as ET
        root = ET.fromstring(fetch_resp.content)
        articles = []
        for article in root.findall(".//PubmedArticle"):
            pmid_el = article.find(".//PMID")
            title_el = article.find(".//ArticleTitle")
            abstract_els = article.findall(".//AbstractText")
            year_el = article.find(".//PubDate/Year")
            journal_el = article.find(".//Journal/Title")

            pmid = pmid_el.text if pmid_el is not None else ""
            title = title_el.text if title_el is not None else ""
            abstract = " ".join(
                (el.text or "") for el in abstract_els
            ).strip()
            year = year_el.text if year_el is not None else ""
            journal = journal_el.text if journal_el is not None else ""

            if abstract:
                articles.append({
                    "pmid": pmid, "title": title,
                    "abstract": abstract, "journal": journal, "year": year,
                })

        logger.info(f"N9: PubMed returned {len(articles)} abstracts for query: {query[:60]}")
        return articles

    except Exception as e:
        logger.warning(f"N9: PubMed fetch failed: {e}")
        return []


def _build_collection_name(run_id: str, protein_id: str) -> str:
    h = hashlib.md5(f"{run_id}_{protein_id}".encode()).hexdigest()[:8]
    return f"td_{h}"


def _index_abstracts_qdrant(
    client, collection_name: str, abstracts: List[Dict], embedder
) -> bool:
    """Index abstracts into Qdrant collection. Returns True on success."""
    try:
        from qdrant_client.models import Distance, VectorParams, PointStruct
        import numpy as np

        texts = [f"{a['title']} {a['abstract']}" for a in abstracts]
        vectors = embedder.encode(texts, show_progress_bar=False).tolist()

        client.recreate_collection(
            collection_name=collection_name,
            vectors_config=VectorParams(size=len(vectors[0]), distance=Distance.COSINE),
        )
        client.upsert(
            collection_name=collection_name,
            points=[
                PointStruct(
                    id=i,
                    vector=vectors[i],
                    payload=abstracts[i],
                )
                for i in range(len(abstracts))
            ],
        )
        return True
    except Exception as e:
        logger.warning(f"N9: Qdrant index failed: {e}")
        return False


def _search_evidence(
    client, collection_name: str, query_text: str, embedder, top_k: int = 3
) -> List[Dict]:
    """Search indexed abstracts for evidence about a specific epitope/query."""
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
        logger.warning(f"N9: Qdrant search failed: {e}")
        return []


def _detect_failure_signals(abstracts: List[Dict]) -> List[str]:
    """
    Scan abstracts for known failure signals relevant to vaccine candidates.
    Returns list of human-readable failure signal strings.
    """
    failure_keywords = {
        "immune evasion":       "immune evasion reported",
        "highly variable":      "high sequence variability escape risk",
        "polymorphic":          "polymorphic region population coverage may vary",
        "cross-reactive":       "cross-reactivity with host proteins reported",
        "poor immunogen":       "poor immunogenicity reported in human trials",
        "no t cell response":   "absence of T-cell response in human subjects",
        "failed phase":         "clinical trial failure reported",
        "toxic":                "toxicity concerns in literature",
        "allergen":             "allergenicity signal in literature",
        "secreted":             "protein is secreted may not be surface-accessible",
    }

    found = []
    for abstract in abstracts:
        text = (abstract.get("abstract", "") + " " + abstract.get("title", "")).lower()
        for keyword, signal in failure_keywords.items():
            if keyword in text and signal not in found:
                found.append(signal)

    return found


def _synthesize_with_claude(
    protein_name: str, abstracts: List[Dict], failure_signals: List[str]
) -> Optional[str]:
    """
    Use Claude to synthesize literature evidence into a one-paragraph summary.
    Only called if ANTHROPIC_API_KEY is set. Non-blocking if unavailable.
    """
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
(2) any known concerns or failure signals, (3) your assessment of the evidence quality.
Be precise and scientific. No hedging language.

{failure_text}

Literature:
{abstract_text}"""

        resp = requests.post(
            "https://api.anthropic.com/v1/messages",
            headers={
                "x-api-key": api_key,
                "anthropic-version": "2023-06-01",
                "content-type": "application/json",
            },
            json={
                "model": "claude-sonnet-4-6",
                "max_tokens": 300,
                "messages": [{"role": "user", "content": prompt}],
            },
            timeout=20,
        )
        resp.raise_for_status()
        return resp.json()["content"][0]["text"].strip()
    except Exception as e:
        logger.warning(f"N9: Claude synthesis failed: {e}")
        return None


class LiteratureAgent:
    """
    Literature Agent
    Searches PubMed for prior evidence on each candidate protein and its predicted epitopes.
    Uses Qdrant for semantic search over abstracts.
    """

    def __init__(self):
        self.stage_name = "literature_search"
        self._qdrant    = None
        self._embedder  = None
        self._use_qdrant = True

    def _init_clients(self):
        if self._qdrant is None:
            self._qdrant = _get_qdrant_client()
            if self._qdrant is None:
                self._qdrant = _get_chroma_client()
                self._use_qdrant = False
        if self._embedder is None:
            self._embedder = _get_embedder()

    def run(
        self,
        candidates: List[CandidateProtein],
        run_id: str = "unknown",
    ) -> List[CandidateProtein]:
        logger.info("N9: Starting literature search")
        self._init_clients()

        active = [c for c in candidates if c.status.value == "active"]
        logger.info(f"   {len(active)} candidates")

        for i, candidate in enumerate(active):
            logger.info(f"   [{i+1}/{len(active)}] {candidate.protein_name}")
            start = time.time()

            # Build PubMed query
            query = (
                f'"{candidate.protein_name}"[Title/Abstract] AND '
                f'("T cell"[Title/Abstract] OR "epitope"[Title/Abstract] OR '
                f'"vaccine"[Title/Abstract]) AND '
                f'"Mycobacterium tuberculosis"[Organism] OR '
                f'"{candidate.protein_id}"[Title/Abstract]'
            )

            abstracts = _fetch_pubmed_abstracts(query)

            # Additional query for specific epitope sequences
            high_conf_epitopes = [
                ep for ep in candidate.ctl_epitopes
                if ep.confidence_tier == ConfidenceTier.HIGH
            ][:3]

            epitope_abstracts = []
            for ep in high_conf_epitopes:
                ep_query = f'"{ep.sequence}"[Title/Abstract] AND "T cell"[Title/Abstract]'
                ep_abs = _fetch_pubmed_abstracts(ep_query, max_results=5)
                epitope_abstracts.extend(ep_abs)

            all_abstracts = abstracts + epitope_abstracts
            # Deduplicate by PMID
            seen_pmids = set()
            unique_abstracts = []
            for a in all_abstracts:
                if a["pmid"] not in seen_pmids:
                    seen_pmids.add(a["pmid"])
                    unique_abstracts.append(a)

            logger.info(f"      {len(unique_abstracts)} unique abstracts")

            # Detect failure signals
            failure_signals = _detect_failure_signals(unique_abstracts)

            # Check prior validation
            prior_validated = any(
                ep.sequence in (a.get("abstract", "") + a.get("title", ""))
                for ep in high_conf_epitopes
                for a in unique_abstracts
            )

            # Index in Qdrant for semantic search
            collection_name = _build_collection_name(run_id, candidate.protein_id)
            if self._qdrant and self._embedder and unique_abstracts and self._use_qdrant:
                indexed = _index_abstracts_qdrant(
                    self._qdrant, collection_name, unique_abstracts, self._embedder
                )
                if indexed:
                    logger.info(f"      Indexed {len(unique_abstracts)} abstracts in Qdrant")

            # Claude synthesis
            literature_summary = _synthesize_with_claude(
                candidate.protein_name, unique_abstracts, failure_signals
            )

            elapsed = time.time() - start

            candidate.add_decision(
                stage=self.stage_name,
                decision="literature_searched",
                reasoning=(
                    f"PubMed search for '{candidate.protein_name}' returned "
                    f"{len(unique_abstracts)} abstracts. "
                    f"Prior experimental validation found: {prior_validated}. "
                    f"Failure signals: {', '.join(failure_signals) if failure_signals else 'none detected'}. "
                    f"Evidence PMIDs: {', '.join(a['pmid'] for a in unique_abstracts[:5])}. "
                    f"Search time: {elapsed:.1f}s."
                    + (f" Literature summary: {literature_summary}" if literature_summary else "")
                ),
                pubmed_hits=len(unique_abstracts),
                prior_validated=prior_validated,
                evidence_pmids=[a["pmid"] for a in unique_abstracts[:10]],
                failure_signals=failure_signals,
                literature_summary=literature_summary,
                search_time_s=round(elapsed, 1),
            )

            logger.info(
                f"      prior_validated={prior_validated} | "
                f"failure_signals={len(failure_signals)} | "
                f"pmids={len(unique_abstracts)} | {elapsed:.1f}s"
            )

        logger.info("N9: Literature search complete")
        return candidates

    def get_status(self) -> Dict[str, Any]:
        return {
            "qdrant_available": self._qdrant is not None,
            "embedder_available": self._embedder is not None,
            "qdrant_url": QDRANT_URL or "in-memory",
            "pubmed_api_key": bool(PUBMED_API_KEY),
        }


literature_agent = LiteratureAgent()