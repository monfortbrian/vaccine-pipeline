"""
IEDB CLIENT - T-CELL EPITOPE PREDICTION
License-free alternative to NetMHCpan using IEDB hosted APIs.

IEDB MHC-I columns: allele seq_num start end length peptide core icore score percentile_rank
IEDB MHC-II columns: allele seq_num start end length core_peptide peptide method percentile_rank

The 'score' is EL presentation score (0-1), NOT IC50.
The 'percentile_rank' is what matters:
  rank < 0.5 = strong binder
  rank < 2.0 = weak binder
"""

import requests
import time
import logging
from typing import List, Dict, Any, Optional

logger = logging.getLogger(__name__)


class IEDBClient:
    def __init__(self):
        self.base_url = "http://tools-cluster-interface.iedb.org/tools_api"
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "Kozi-Pipeline/2.0"})
        self.max_retries = 3

        self.mhc_i_alleles = [
            "HLA-A*01:01", "HLA-A*02:01", "HLA-A*03:01", "HLA-A*11:01", "HLA-A*24:02",
            "HLA-A*30:01", "HLA-A*68:01", "HLA-B*07:02", "HLA-B*08:01", "HLA-B*15:01",
            "HLA-B*35:01", "HLA-B*40:01", "HLA-B*44:02", "HLA-B*51:01", "HLA-B*53:01"
        ]

        self.mhc_ii_alleles = [
            "DRB1*01:01", "DRB1*03:01", "DRB1*04:01", "DRB1*07:01",
            "DRB1*11:01", "DRB1*13:01", "DRB1*15:01"
        ]

    def predict_mhc_i_binding(self, sequence: str, alleles: Optional[List[str]] = None) -> List[Dict[str, Any]]:
        """Predict CTL epitopes. Chunks alleles to avoid IEDB timeouts."""
        if alleles is None:
            alleles = self.mhc_i_alleles
        all_epitopes = []
        chunk_size = 3
        for length in [9, 10]:
            for i in range(0, len(alleles), chunk_size):
                chunk = alleles[i:i+chunk_size]
                try:
                    epitopes = self._call_with_retry(
                        f"{self.base_url}/mhci/",
                        {
                            'method': 'netmhcpan_el',
                            'sequence_text': f">query\n{sequence}",
                            'allele': ','.join(chunk),
                            'length': ','.join([str(length)] * len(chunk))
                        },
                        'CTL', length
                    )
                    all_epitopes.extend(epitopes)
                    time.sleep(1.5)
                except Exception as e:
                    logger.warning(f"CTL prediction failed for {chunk}: {e}")
        all_epitopes.sort(key=lambda x: x.get('ic50_nm', 50000))
        logger.info(f"Predicted {len(all_epitopes)} CTL epitopes")
        return all_epitopes

    def predict_mhc_ii_binding(self, sequence: str, alleles: Optional[List[str]] = None) -> List[Dict[str, Any]]:
        """Predict HTL epitopes. Chunks alleles to avoid IEDB timeouts."""
        if alleles is None:
            alleles = self.mhc_ii_alleles
        all_epitopes = []
        chunk_size = 3
        for length in [15]:
            for i in range(0, len(alleles), chunk_size):
                chunk = alleles[i:i+chunk_size]
                try:
                    epitopes = self._call_with_retry(
                        f"{self.base_url}/mhcii/",
                        {
                            'method': 'netmhciipan',
                            'sequence_text': f">query\n{sequence}",
                            'allele': ','.join(chunk),
                            'length': ','.join([str(length)] * len(chunk))
                        },
                        'HTL', length
                    )
                    all_epitopes.extend(epitopes)
                    time.sleep(1.5)
                except Exception as e:
                    logger.warning(f"HTL prediction failed for {chunk}: {e}")
        all_epitopes.sort(key=lambda x: x.get('ic50_nm', 50000))
        logger.info(f"Predicted {len(all_epitopes)} HTL epitopes")
        return all_epitopes

    def _call_with_retry(self, url: str, data: dict, epitope_type: str, length: int) -> List[Dict[str, Any]]:
        """POST to IEDB with exponential backoff retry on failure."""
        for attempt in range(self.max_retries):
            try:
                response = self.session.post(url, data=data, timeout=120)
                response.raise_for_status()

                # IEDB sometimes returns HTML error pages with 200 status
                if "<html" in response.text.lower()[:200]:
                    raise ValueError("IEDB returned HTML error page instead of TSV")

                return self._parse_mhc_results(response.text, epitope_type, length)

            except Exception as e:
                wait = (2 ** attempt) * 2  # 2s, 4s, 8s
                if attempt < self.max_retries - 1:
                    logger.debug(f"  IEDB retry {attempt+1}/{self.max_retries}: {e}")
                    time.sleep(wait)
                else:
                    raise

    def _parse_mhc_results(self, response_text: str, epitope_type: str, length: int) -> List[Dict[str, Any]]:
        """Parse IEDB response using header-based column detection."""
        epitopes = []
        try:
            lines = response_text.strip().split('\n')
            if len(lines) < 2:
                return epitopes

            header = lines[0].split('\t')
            col = {}
            for idx, name in enumerate(header):
                col[name.strip().lower()] = idx

            allele_idx = col.get('allele', 0)
            peptide_idx = col.get('peptide', 5)
            start_idx = col.get('start', 2)
            end_idx = col.get('end', 3)
            score_idx = col.get('score')
            rank_idx = col.get('percentile_rank')

            if rank_idx is None:
                for alt in ['rank', 'percentile', '%rank']:
                    if alt in col:
                        rank_idx = col[alt]
                        break

            for line in lines[1:]:
                if not line.strip():
                    continue
                parts = line.split('\t')
                if len(parts) < 6:
                    continue
                try:
                    allele = parts[allele_idx].strip()
                    peptide = parts[peptide_idx].strip()
                    if not peptide or len(peptide) < 5:
                        continue

                    start_pos = int(parts[start_idx].strip()) if start_idx < len(parts) else 0
                    end_pos = int(parts[end_idx].strip()) if end_idx < len(parts) else 0

                    score = None
                    if score_idx is not None and score_idx < len(parts):
                        try:
                            score = float(parts[score_idx].strip())
                        except ValueError:
                            pass

                    rank = None
                    if rank_idx is not None and rank_idx < len(parts):
                        try:
                            rank = float(parts[rank_idx].strip())
                        except ValueError:
                            pass

                    if rank is not None:
                        if rank <= 0.5:
                            approx_ic50 = rank * 100
                        elif rank <= 2.0:
                            approx_ic50 = 50 + (rank - 0.5) * 300
                        elif rank <= 10:
                            approx_ic50 = 500 + (rank - 2.0) * 562.5
                        else:
                            approx_ic50 = 5000 + rank * 100
                    elif score is not None:
                        approx_ic50 = (1 - score) * 5000 if score > 0 else 50000
                    else:
                        approx_ic50 = 50000

                    epitopes.append({
                        'sequence': peptide,
                        'allele': allele,
                        'ic50_nm': round(approx_ic50, 1),
                        'percentile_rank': rank,
                        'el_score': score,
                        'length': len(peptide),
                        'epitope_type': epitope_type,
                        'prediction_method': f'IEDB_{epitope_type}',
                        'strong_binder': rank is not None and rank < 0.5,
                        'weak_binder': rank is not None and rank < 2.0,
                        'start_position': start_pos,
                        'end_position': end_pos,
                    })
                except (ValueError, IndexError):
                    continue

            max_ic50 = 5000 if epitope_type == 'CTL' else 10000
            epitopes = [e for e in epitopes if e['ic50_nm'] < max_ic50]
        except Exception as e:
            logger.error(f"Failed to parse IEDB {epitope_type} results: {e}")
        return epitopes

    def test_connection(self) -> bool:
        try:
            response = self.session.post(
                f"{self.base_url}/mhci/",
                data={
                    'method': 'netmhcpan_el',
                    'sequence_text': ">test\nMKLRLFCLAMLMACAQILNGS",
                    'allele': 'HLA-A*02:01',
                    'length': '9'
                },
                timeout=30
            )
            return response.status_code == 200 and len(response.text) > 50
        except Exception:
            return False


iedb = IEDBClient()