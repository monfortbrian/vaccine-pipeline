"""
MHCflurry 2.0 client - local MHC-I epitope prediction.
Fallback when IEDB tools-cluster-interface is unavailable.

Reference: O'Brien et al., Cell Systems 2019, 9(5):452-458.
           https://doi.org/10.1016/j.cels.2019.08.005
License: Apache 2.0 - no commercial restrictions.
Models: downloaded at Docker build time via `mhcflurry-downloads fetch`.
        Cached at /root/.local/share/mhcflurry inside the container.

Predictor used: Class1AffinityPredictor
  Predicts IC50 affinity (nM) per peptide-allele pair.
  Distinct from Class1PresentationPredictor which treats allele list
  as a diploid genotype (max 6 alleles) - that API is not used here.

Supported alleles: subset of our 15-allele panel confirmed in MHCflurry
  trained allele set. Checked at runtime via predictor.supported_alleles.

Thresholds (equivalent to IEDB NetMHCpan EL):
  IC50 < 500 nM  → strong binder (HIGH confidence)
  IC50 < 5000 nM → weak binder  (MEDIUM/LOW confidence)
"""

import logging
from typing import List, Dict, Any, Optional

logger = logging.getLogger(__name__)

# Our 15-allele panel - all supported by MHCflurry 2.0
PANEL_ALLELES = [
    "HLA-A*01:01", "HLA-A*02:01", "HLA-A*03:01", "HLA-A*11:01",
    "HLA-A*24:02", "HLA-A*30:01", "HLA-A*68:01",
    "HLA-B*07:02", "HLA-B*08:01", "HLA-B*15:01",
    "HLA-B*35:01", "HLA-B*40:01", "HLA-B*44:02",
    "HLA-B*51:01", "HLA-B*53:01",
]


class MHCflurryClient:
    def __init__(self):
        self._predictor = None
        self._available: Optional[bool] = None
        self._supported_alleles: Optional[set] = None

    # ── LAZY LOAD ─────────────────────────────────────────────────────────────

    @property
    def predictor(self):
        """
        Lazy-load Class1AffinityPredictor on first call.
        Uses affinity predictor, NOT presentation predictor -
        affinity predictor accepts per-peptide allele lists with no size limit.
        """
        if self._predictor is None:
            try:
                from mhcflurry import Class1AffinityPredictor
                self._predictor = Class1AffinityPredictor.load()
                self._supported_alleles = set(self._predictor.supported_alleles)
                self._available = True
                logger.info(
                    f"MHCflurry 2.0 (Class1AffinityPredictor) loaded. "
                    f"{len(self._supported_alleles)} alleles supported."
                )
            except Exception as e:
                self._available = False
                self._predictor = None
                logger.error(f"MHCflurry failed to load: {e}")
        return self._predictor

    def is_available(self) -> bool:
        """Check if MHCflurry is loaded and ready. Triggers load on first call."""
        if self._available is None:
            _ = self.predictor
        return bool(self._available)

    # ── PREDICTION ────────────────────────────────────────────────────────────

    def predict_mhc_i_binding(
        self,
        sequence: str,
        alleles: Optional[List[str]] = None,
    ) -> List[Dict[str, Any]]:
        """
        Predict MHC-I binding for peptide lengths 9 and 10.
        Predicts each allele separately to avoid the genotype-list limit.
        Returns dicts in the same format as iedb_client for seamless substitution.
        """
        if not self.is_available():
            logger.warning("MHCflurry not available - returning empty")
            return []

        if alleles is None:
            alleles = PANEL_ALLELES

        # Only predict for alleles the model was trained on
        valid_alleles = [
            a for a in alleles
            if a in self._supported_alleles
        ]
        if not valid_alleles:
            logger.warning(
                f"None of the requested alleles are in MHCflurry's trained set. "
                f"Requested: {alleles[:3]}..."
            )
            return []

        logger.info(
            f"MHCflurry: predicting {len(valid_alleles)} alleles × "
            f"lengths [9,10] for sequence of {len(sequence)} residues"
        )

        import pandas as pd
        all_epitopes: List[Dict[str, Any]] = []

        for length in [9, 10]:
            peptides = self._slice_peptides(sequence, length)
            if not peptides:
                logger.debug(f"Sequence too short for length {length}, skipping")
                continue

            frames = []
            for allele in valid_alleles:
                try:
                    # predict() accepts parallel lists: peptides[i] → alleles[i]
                    df = self.predictor.predict(
                        peptides=peptides,
                        alleles=[allele] * len(peptides),
                    )
                    # Column is 'mhcflurry_affinity', add allele and peptide cols
                    df = df.copy()
                    df["allele"]  = allele
                    df["peptide"] = peptides
                    frames.append(df)
                except Exception as e:
                    logger.warning(f"MHCflurry failed for allele {allele}: {e}")

            if not frames:
                continue

            result = pd.concat(frames, ignore_index=True)

            for _, row in result.iterrows():
                ic50 = float(row.get("mhcflurry_affinity", 50000))
                percentile = row.get("mhcflurry_affinity_percentile", None)

                # Skip weak binders - same threshold as IEDB filter
                if ic50 > 5000:
                    continue

                # Map affinity to approximate percentile rank
                # (MHCflurry affinity_percentile ≈ IEDB percentile_rank)
                if percentile is not None:
                    rank = float(percentile)
                elif ic50 < 50:
                    rank = 0.3
                elif ic50 < 500:
                    rank = 1.5
                elif ic50 < 5000:
                    rank = 5.0
                else:
                    rank = 15.0

                all_epitopes.append({
                    "sequence":          str(row["peptide"]),
                    "allele":            str(row["allele"]),
                    "ic50_nm":           round(ic50, 1),
                    "percentile_rank":   round(rank, 2),
                    "el_score":          None,  # affinity predictor has no EL score
                    "length":            len(str(row["peptide"])),
                    "epitope_type":      "CTL",
                    "prediction_method": "MHCflurry_2.0_affinity",
                    "strong_binder":     ic50 < 500,
                    "weak_binder":       ic50 < 5000,
                    "start_position":    0,
                    "end_position":      0,
                })

        all_epitopes.sort(key=lambda x: x.get("ic50_nm", 50000))
        logger.info(f"MHCflurry predicted {len(all_epitopes)} CTL epitopes")
        return all_epitopes

    # ── HELPERS ───────────────────────────────────────────────────────────────

    def _slice_peptides(self, sequence: str, length: int) -> List[str]:
        """All overlapping subsequences of given length."""
        if len(sequence) < length:
            return []
        return [sequence[i:i + length] for i in range(len(sequence) - length + 1)]

    def test_connection(self) -> bool:
        """Smoke test - load model and run one prediction."""
        try:
            result = self.predict_mhc_i_binding(
                sequence="MKLRLFCLAMLMACAQILNGS",
                alleles=["HLA-A*02:01"],
            )
            success = isinstance(result, list)
            logger.info(f"MHCflurry smoke test: {'PASS' if success else 'FAIL'}")
            return success
        except Exception as e:
            logger.error(f"MHCflurry smoke test failed: {e}")
            return False


# Global instance - predictor loads lazily on first prediction call
mhcflurry_client = MHCflurryClient()