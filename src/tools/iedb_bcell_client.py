"""
IEDB B-CELL CLIENT - ANTIBODY EPITOPE PREDICTION
Uses IEDB BepiPred 2.0 - license-free alternative to BepiPred 3.0.

Response formats:
  BepiPred: Position Residue Score Assignment  (score = col 2)
  Parker:   Position Residue Start End Peptide Score  (score = last col)
  Emini:    Position Residue Start End Peptide Score  (score = last col)
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

    def predict_linear_epitopes(self, sequence, method='bepipred2'):
        try:
            scores = self._call_iedb_bcell_api(sequence, method)
            scores = self._pad_scores(scores, len(sequence), method)
            epitopes = self._extract_epitope_regions(sequence, scores, method)
            logger.info(f"B-cell prediction: {len(epitopes)} epitopes found using {method}")
            return epitopes
        except Exception as e:
            logger.error(f"B-cell prediction failed: {e}")
            return []

    def predict_consensus_epitopes(self, sequence):
        all_predictions = {}
        seq_len = len(sequence)
        for method_key in ['bepipred2', 'parker', 'emini']:
            try:
                predictions = self._call_iedb_bcell_api(sequence, method_key)
                if predictions and len(predictions) > 0:
                    predictions = self._pad_scores(predictions, seq_len, method_key)
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

    def _pad_scores(self, scores, target_len, method):
        if len(scores) == target_len:
            return scores
        if len(scores) == 0:
            return scores
        gap = target_len - len(scores)
        if gap <= 0:
            return scores[:target_len]
        pad_left = gap // 2
        pad_right = gap - pad_left
        padded = ([scores[0]] * pad_left) + scores + ([scores[-1]] * pad_right)
        return padded[:target_len]

    def _call_iedb_bcell_api(self, sequence, method):
        method_name = self.methods.get(method, 'Bepipred-2.0')
        data = {
            'method': method_name,
            'sequence_text': f">query\n{sequence}",
        }
        if method != 'bepipred2':
            data['window_size'] = '7'
        for attempt in range(3):
            try:
                response = self.session.post(
                    f"{self.base_url}/bcell/",
                    data=data,
                    timeout=120
                )
                response.raise_for_status()
                return self._parse_iedb_bcell_response(response.text, method)
            except Exception as e:
                if attempt < 2:
                    time.sleep((attempt + 1) * 2)
                else:
                    raise

    def _parse_iedb_bcell_response(self, response_text, method):
        scores = []
        try:
            lines = response_text.strip().split('\n')
            for line in lines:
                line = line.strip()
                if not line:
                    continue
                if line.lower().startswith('position') or line.lower().startswith('pos'):
                    continue
                parts = line.split()
                if len(parts) < 3:
                    continue
                try:
                    position = int(parts[0])
                    residue = parts[1]
                    if not (len(residue) == 1 and residue.isalpha()):
                        continue
                    # Find last numeric value in the row - that is always the score
                    # BepiPred: Position Residue Score Assignment -> score is parts[2]
                    # Parker/Emini: Position Residue Start End Peptide Score -> score is parts[-1]
                    score = None
                    for p in reversed(parts[2:]):
                        try:
                            val = float(p)
                            score = val
                            break
                        except ValueError:
                            continue
                    if score is not None:
                        scores.append(score)
                except (ValueError, IndexError):
                    continue
            logger.debug(f"Parsed {len(scores)} B-cell scores from IEDB {method}")
        except Exception as e:
            logger.error(f"Error parsing IEDB B-cell response: {e}")
        return scores

    def _extract_epitope_regions(self, sequence, scores, method):
        if not scores or len(scores) != len(sequence):
            logger.warning(f"Score length mismatch: {len(scores)} vs {len(sequence)}")
            return []
        epitopes = []
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
                if 6 <= length <= 25:
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

    def _find_consensus_regions(self, sequence, all_predictions):
        seq_len = len(sequence)
        consensus_scores = [0.0] * seq_len
        method_count = 0
        for method, scores in all_predictions.items():
            if len(scores) == seq_len:
                if method == 'bepipred2':
                    norm_scores = scores
                else:
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

    def _score_confidence(self, score, method):
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

    def test_connection(self):
        try:
            result = self.predict_linear_epitopes("MKLRLFCLAMLMACAQILNGS")
            return isinstance(result, list)
        except Exception:
            return False


iedb_bcell = IEDBBCellClient()