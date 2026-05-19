"""
MHCflurry 2.0 client - local MHC-I epitope prediction.
Used as primary when IEDB legacy API is unavailable.

Reference: O'Brien et al., Cell Systems 2019, 9(5):452-458.
License: Apache 2.0 - no commercial restrictions.
Models: downloaded at Docker build time via mhcflurry-downloads fetch.

Supported alleles: 140+ HLA-A/B/C alleles.
Output: affinity (IC50 nM) + presentation score (0-1).
Threshold: IC50 < 500 nM = strong binder (equivalent to IEDB rank < 2.0)
"""

import logging
from typing import List, Dict, Any, Optional

logger = logging.getLogger(__name__)

# Alleles supported by MHCflurry that overlap with our panel
SUPPORTED_ALLELES = {
    "HLA-A*01:01", "HLA-A*02:01", "HLA-A*03:01", "HLA-A*11:01",
    "HLA-A*24:02", "HLA-A*30:01", "HLA-A*68:01",
    "HLA-B*07:02", "HLA-B*08:01", "HLA-B*15:01",
    "HLA-B*35:01", "HLA-B*40:01", "HLA-B*44:02",
    "HLA-B*51:01", "HLA-B*53:01",
}


class MHCflurryClient:
    def __init__(self):
        self._predictor = None
        self._available = None

    @property
    def predictor(self):
        """Lazy load - only import when first needed."""
        if self._predictor is None:
            try:
                from mhcflurry import Class1PresentationPredictor
                self._predictor = Class1PresentationPredictor.load()
                self._available = True
                logger.info("MHCflurry 2.0 loaded successfully")
            except Exception as e:
                self._available = False
                logger.error(f"MHCflurry failed to load: {e}")
        return self._predictor

    def is_available(self) -> bool:
        """Check if MHCflurry is loaded and ready."""
        if self._available is None:
            _ = self.predictor
        return self._available or False

    def predict_mhc_i_binding(
        self,
        sequence: str,
        alleles: Optional[List[str]] = None,
    ) -> List[Dict[str, Any]]:
        """
        Predict MHC-I binding for all peptide lengths 9-10.
        Returns list of dicts compatible with iedb_client output format.
        """
        if not self.is_available():
            logger.warning("MHCflurry not available - returning empty")
            return []

        if alleles is None:
            alleles = list(SUPPORTED_ALLELES)

        # Filter to alleles MHCflurry supports
        valid_alleles = [a for a in alleles if a in SUPPORTED_ALLELES]
        if not valid_alleles:
            logger.warning("No valid MHCflurry alleles in panel")
            return []

        all_epitopes = []

        for length in [9, 10]:
            peptides = self._slice_peptides(sequence, length)
            if not peptides:
                continue

            try:
                result = self.predictor.predict(
                    peptides=peptides,
                    alleles=valid_alleles,
                    include_affinity_percentile=True,
                )

                for _, row in result.iterrows():
                    ic50 = row.get("affinity", 50000)
                    presentation_score = row.get("presentation_score", 0)
                    percentile = row.get("affinity_percentile", None)

                    # Skip weak binders
                    if ic50 > 5000:
                        continue

                    # Approximate percentile rank from affinity
                    if percentile is not None:
                        rank = percentile
                    elif ic50 < 50:
                        rank = 0.3
                    elif ic50 < 500:
                        rank = 1.5
                    elif ic50 < 5000:
                        rank = 5.0
                    else:
                        rank = 15.0

                    all_epitopes.append({
                        "sequence":          row["peptide"],
                        "allele":            row["allele"],
                        "ic50_nm":           round(float(ic50), 1),
                        "percentile_rank":   round(float(rank), 2) if rank else None,
                        "el_score":          round(float(presentation_score), 4),
                        "length":            len(row["peptide"]),
                        "epitope_type":      "CTL",
                        "prediction_method": "MHCflurry_2.0",
                        "strong_binder":     ic50 < 500,
                        "weak_binder":       ic50 < 5000,
                        "start_position":    0,
                        "end_position":      0,
                    })

            except Exception as e:
                logger.warning(f"MHCflurry prediction failed for length {length}: {e}")

        all_epitopes.sort(key=lambda x: x.get("ic50_nm", 50000))
        logger.info(f"MHCflurry predicted {len(all_epitopes)} CTL epitopes")
        return all_epitopes

    def _slice_peptides(self, sequence: str, length: int) -> List[str]:
        """Slice all overlapping peptides of given length from sequence."""
        if len(sequence) < length:
            return []
        return [sequence[i:i+length] for i in range(len(sequence) - length + 1)]

    def test_connection(self) -> bool:
        """Verify MHCflurry is loaded and can make a prediction."""
        try:
            result = self.predict_mhc_i_binding(
                "MKLRLFCLAMLMACAQILNGS",
                alleles=["HLA-A*02:01"],
            )
            return isinstance(result, list)
        except Exception:
            return False


# Global instance
mhcflurry_client = MHCflurryClient()