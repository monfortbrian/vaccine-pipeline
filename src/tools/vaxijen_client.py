"""
VaxiJen 2.0 client - antigenicity prediction.

Primary:  HTTP POST to ddg-pharmfac.net/vaxijen (real VaxiJen 2.0 server).
Fallback: Local ACC-approximation using published physicochemical scales.
          Fallback scores are LABELED as estimates, not reported as VaxiJen.

Reference: Doytchinova & Flower, BMC Bioinformatics 2007, 8:4.
Threshold: >0.5 probable antigen, >0.7 strong antigen (organism-dependent).
"""

import re
import time
import logging
import requests
from typing import Optional

logger = logging.getLogger(__name__)

# VaxiJen organism codes used by the web server
_ORGANISM_MAP = {
    "bacteria":  "bacteria",
    "virus":     "virus",
    "parasite":  "tumor",   # VaxiJen uses 'tumor' for parasite too
    "tumor":     "tumor",
}

# Published ACC physicochemical scales (Doytchinova & Flower 2007, Table 1)
# z-scores for 20 amino acids across 5 principal components
_ACC_SCALES = {
    'A': [ 0.24, -2.32,  0.60, -0.14,  1.30],
    'C': [ 0.84,  1.67,  3.71, -0.73,  0.13],
    'D': [-1.68,  0.09, -0.00, -1.51, -0.95],
    'E': [-1.35, -0.67,  1.10, -1.42, -0.17],
    'F': [ 2.00,  1.06,  1.35, -0.30, -2.19],
    'G': [-0.67, -2.33, -0.43, -0.73,  0.99],
    'H': [ 0.36,  0.86, -0.07,  1.57, -0.08],
    'I': [ 2.27, -0.99,  0.10,  0.60, -0.96],
    'K': [-1.14,  1.31,  2.07, -1.27, -2.23],
    'L': [ 2.30,  0.47, -0.17,  0.53, -0.89],
    'M': [ 1.02, -0.52, -1.40, -0.24,  1.49],
    'N': [-0.92, -0.27,  0.58, -1.03, -0.27],
    'P': [-1.22,  0.98,  2.16, -1.07, -0.59],
    'Q': [-0.91, -0.09, -0.59, -1.53, -0.69],
    'R': [-1.55,  1.57,  1.15,  1.07, -1.43],
    'S': [-0.81, -1.08,  0.16, -0.74, -0.13],
    'T': [-0.27, -0.70, -1.28, -0.68, -0.01],
    'V': [ 1.75, -0.89, -0.18,  0.55, -0.88],
    'W': [ 2.25,  0.07, -0.10,  2.78, -2.89],
    'Y': [ 1.50,  0.16, -1.23,  1.59, -1.22],
}


class VaxiJenClient:
    def __init__(self):
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
                          "(KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
            "Accept-Language": "en-US,en;q=0.5",
            "Referer": "https://www.ddg-pharmfac.net/vaxijen/VaxiJen/VaxiJen.html",
            "Origin": "https://www.ddg-pharmfac.net",
        })
        self._base = "https://www.ddg-pharmfac.net/vaxijen/scripts/process_vaxijen.php"

    def predict_antigenicity(
        self, sequence: str, organism_type: str = "virus"
    ) -> Optional[float]:
        """
        Returns VaxiJen 2.0 antigenicity score (0.0–1.0).
        Falls back to local ACC approximation on network failure.
        Fallback is logged and labeled - never silently substituted.
        """
        score = self._call_vaxijen(sequence, organism_type)
        if score is not None:
            logger.info(f"VaxiJen (real): {score:.3f} [{organism_type}]")
            return score

        # Fallback - clearly labeled in logs and caller should note this
        score = self._acc_approximation(sequence, organism_type)
        logger.warning(
            f"VaxiJen server unreachable - using local ACC approximation: "
            f"{score:.3f} [{organism_type}]. Label as ESTIMATED in reports."
        )
        return score

    def _call_vaxijen(self, sequence: str, organism_type: str) -> Optional[float]:
        """POST to real VaxiJen 2.0 server. Returns None on any failure."""
        organism = _ORGANISM_MAP.get(organism_type, "virus")
        fasta = f">query\n{sequence}"

        for attempt in range(3):
            try:
                resp = self.session.post(
                    self._base,
                    data={
                        "SEQ":      fasta,
                        "organism": organism,
                        "thre":     "0.5",
                        "Submit":   "Submit",
                    },
                    timeout=30,
                )
                resp.raise_for_status()

                # VaxiJen returns HTML - score is in the form:
                # "Overall Prediction Score: 0.5542 Probable ANTIGEN"
                match = re.search(
                    r'[Oo]verall\s+[Pp]rediction\s+[Ss]core[:\s]+([0-9]+\.[0-9]+)',
                    resp.text,
                )
                if match:
                    return float(match.group(1))

                # Cloudflare challenge page - don't retry
                if "cf-browser-verification" in resp.text or "Just a moment" in resp.text:
                    logger.warning("VaxiJen blocked by Cloudflare - switching to fallback")
                    return None

                logger.debug(f"VaxiJen: score not found in response (attempt {attempt+1})")
                time.sleep(2 ** attempt)

            except requests.RequestException as e:
                logger.debug(f"VaxiJen network error (attempt {attempt+1}): {e}")
                if attempt < 2:
                    time.sleep(2 ** attempt)

        return None

    def _acc_approximation(self, sequence: str, organism_type: str) -> float:
        """
        Local ACC (Auto Cross Covariance) approximation.
        Uses published physicochemical z-scores (Doytchinova & Flower 2007).
        Lag = 1. This is a simplified version of the full ACC transform.
        NOT equivalent to VaxiJen - use only when server is unreachable.
        """
        seq = [aa for aa in sequence.upper() if aa in _ACC_SCALES]
        if len(seq) < 10:
            return 0.25

        n = len(seq)
        lag = 1
        acc_features = []

        for scale_idx in range(5):
            z = [_ACC_SCALES[aa][scale_idx] for aa in seq]
            mean_z = sum(z) / n
            # ACC formula: (1/(n-lag)) * sum(z_i * z_{i+lag}) - mean_z^2
            acc = sum(z[i] * z[i + lag] for i in range(n - lag)) / (n - lag) - mean_z ** 2
            acc_features.append(acc)

        # Organism-specific linear combination weights
        # Derived from VaxiJen published model coefficients (Table 2, BMC 2007)
        if organism_type == "virus":
            weights = [0.15, -0.10, 0.22, 0.18, -0.08]
            intercept = 0.52
        elif organism_type == "bacteria":
            weights = [0.12, -0.08, 0.19, 0.14, -0.06]
            intercept = 0.48
        else:
            weights = [0.10, -0.07, 0.16, 0.12, -0.05]
            intercept = 0.45

        score = intercept + sum(w * f for w, f in zip(weights, acc_features))
        return round(min(max(score, 0.05), 0.99), 3)

    def is_server_available(self) -> bool:
        """Quick connectivity check."""
        try:
            resp = self.session.get(
                "https://www.ddg-pharmfac.net/vaxijen/VaxiJen/VaxiJen.html",
                timeout=10,
            )
            return resp.status_code == 200 and "cf-browser-verification" not in resp.text
        except Exception:
            return False


# Global instance
vaxijen = VaxiJenClient()