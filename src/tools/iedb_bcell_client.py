"""
IEDB B-CELL CLIENT - ANTIBODY EPITOPE PREDICTION
Uses IEDB BepiPred 2.0 - license-free alternative to BepiPred 3.0.

IEDB B-cell response format (tab or space separated):
  Position  Residue  Score  Assignment
  0         M        0.291  .
  1         T        0.360  .
  4         Q        0.548  E      <-- E means epitope (above threshold)

Score range: ~0.0 to ~1.0 (higher = more likely epitope)
Default threshold: 0.5 (scores above = predicted epitope)
Assignment: E = epitope, . = non-epitope
"""

import requests
import logging
import time
from typing import List, Dict, Any

logger = logging.getLogger(__name__)


class IEDBBCellClient:
    def __init__(self):
        self.base_url = "http://tools-cluster-interface.iedb.org/tools_api"
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "Kozi-Pipeline/2.0"})

        self.methods = {
            'bepipred2': 'Bepipred-2.0',
            'parker': 'Parker',
            'emini': 'Emini'
        }

    def predict_linear_epitopes(self, sequence: str, method: str = 'bepipred2') -> List[Dict[str, Any]]:
        """Predict linear B-cell epitopes using a single method."""
        try:
            scores = self._call_iedb_bcell_api(sequence, method)
            epitopes = self._extract_epitope_regions(sequence, scores, method)
            logger.info(f"B-cell prediction: {len(epitopes)} epitopes found using {method}")
            return epitopes
        except Exception as e:
            logger.error(f"B-cell prediction failed: {e}")
            return []

    def predict_consensus_epitopes(self, sequence: str) -> List[Dict[str, Any]]:
        """Use multiple methods for consensus prediction."""
        all_predictions = {}

        for method_key in ['bepipred2', 'parker', 'emini']:
            try:
                predictions = self._call_iedb_bcell_api(sequence, method_key)
                if predictions and len(predictions) > 0:
                    all_predictions[method_key] = predictions
                    logger.debug(f"  {method_key}: {len(predictions)} scores")
                time.sleep(1)
            except Exception as e:
                logger.warning(f"Method {method_key} failed: {e}")

        if not all_predictions:
            return []

        consensus_epitopes = self._find_consensus_regions(sequence, all_predictions)
        logger.info(f"Consensus B-cell prediction: {len(consensus_epitopes)} epitopes")
        return consensus_epitopes

    def _call_iedb_bcell_api(self, sequence: str, method: str) -> List[float]:
        """Call IEDB B-cell prediction API."""
        method_name = self.methods.get(method, 'Bepipred-2.0')

        data = {
            'method': method_name,
            'sequence_text': f">query\n{sequence}",
        }

        # Only add window_size for non-BepiPred methods
        if method != 'bepipred2':
            data['window_size'] = '7'

        response = self.session.post(
            f"{self.base_url}/bcell/",
            data=data,
            timeout=120
        )

        response.raise_for_status()
        return self._parse_iedb_bcell_response(response.text, method)

    def _parse_iedb_bcell_response(self, response_text: str, method: str) -> List[float]:
        """
        Parse IEDB B-cell API response.

        Response format (tab or space separated):
          Position  Residue  Score  Assignment
          0         M        0.291  .
          1         T        0.360  .

        The Score column is what we need - one score per residue.
        """
        scores = []

        try:
            lines = response_text.strip().split('\n')

            for line in lines:
                line = line.strip()
                if not line:
                    continue

                # Skip header line
                if line.lower().startswith('position') or line.lower().startswith('pos'):
                    continue

                # Split by any whitespace (tabs or spaces)
                parts = line.split()

                if len(parts) < 3:
                    continue

                try:
                    # Column 0 = Position (int), Column 1 = Residue (char), Column 2 = Score (float)
                    position = int(parts[0])
                    residue = parts[1]
                    score = float(parts[2])

                    # Validate: residue should be a single amino acid letter
                    if len(residue) == 1 and residue.isalpha():
                        scores.append(score)
                except (ValueError, IndexError):
                    continue

            logger.debug(f"Parsed {len(scores)} B-cell scores from IEDB {method}")

        except Exception as e:
            logger.error(f"Error parsing IEDB B-cell response: {e}")

        return scores

    def _extract_epitope_regions(self, sequence: str, scores: List[float], method: str) -> List[Dict[str, Any]]:
        """Extract epitope regions from prediction scores."""
        if not scores or len(scores) != len(sequence):
            logger.warning(f"Score length mismatch: {len(scores)} vs {len(sequence)}")
            return []

        epitopes = []

        # Use 0.5 threshold for BepiPred-2.0 and consensus; higher for other methods
        if method in ('bepipred2', 'consensus'):
            threshold = 0.5
        else:
            threshold = 1.0

        i = 0
        while i < len(scores):
            if scores[i] > threshold:
                start = i

                while i < len(scores) and scores[i] > threshold:
                    i += 1

                end = i
                length = end - start

                if 6 <= length <= 25:  # Reasonable B-cell epitope length
                    epitope_seq = sequence[start:end]
                    avg_score = sum(scores[start:end]) / length

                    epitope = {
                        'sequence': epitope_seq,
                        'start_position': start + 1,
                        'end_position': end,
                        'length': length,
                        'epitope_score': avg_score,
                        'epitope_type': 'linear',
                        'prediction_method': f'IEDB_{method}',
                        'confidence': self._score_confidence(avg_score, method),
                        'raw_scores': scores[start:end]
                    }
                    epitopes.append(epitope)
            else:
                i += 1

        epitopes.sort(key=lambda x: x['epitope_score'], reverse=True)
        return epitopes[:15]

    def _find_consensus_regions(self, sequence: str, all_predictions: Dict[str, List[float]]) -> List[Dict[str, Any]]:
        """Find consensus epitope regions across multiple methods."""
        seq_len = len(sequence)
        consensus_scores = [0.0] * seq_len
        method_count = 0

        for method, scores in all_predictions.items():
            if len(scores) == seq_len:
                # Normalize scores to 0-1 range for fair comparison
                if method == 'bepipred2':
                    # BepiPred scores are already roughly 0-1
                    norm_scores = scores
                else:
                    # Parker/Emini scores can be negative or > 1, normalize
                    min_s = min(scores) if scores else 0
                    max_s = max(scores) if scores else 1
                    range_s = max_s - min_s if max_s != min_s else 1
                    norm_scores = [(s - min_s) / range_s for s in scores]

                for i in range(seq_len):
                    consensus_scores[i] += norm_scores[i]
                method_count += 1
            else:
                logger.warning(f"Skipping {method}: {len(scores)} scores vs {seq_len} residues")

        if method_count > 0:
            consensus_scores = [score / method_count for score in consensus_scores]

        return self._extract_epitope_regions(sequence, consensus_scores, 'consensus')

    def _score_confidence(self, score: float, method: str) -> str:
        """Assign confidence level."""
        if method in ('bepipred2', 'consensus'):
            if score > 0.7:
                return 'high'
            elif score > 0.55:
                return 'medium'
            elif score > 0.5:
                return 'low'
            else:
                return 'uncertain'
        else:
            if score > 2.0:
                return 'high'
            elif score > 1.5:
                return 'medium'
            elif score > 1.0:
                return 'low'
            else:
                return 'uncertain'

    def test_connection(self) -> bool:
        """Test IEDB B-cell API connection."""
        try:
            test_sequence = "MKLRLFCLAMLMACAQILNGS"
            result = self.predict_linear_epitopes(test_sequence)
            success = isinstance(result, list)
            logger.info(f"IEDB B-cell test: {'Success' if success else 'Failed'}")
            return success
        except Exception as e:
            logger.error(f"IEDB B-cell test failed: {e}")
            return False


# Global instance
iedb_bcell = IEDBBCellClient()