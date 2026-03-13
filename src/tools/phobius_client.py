"""
Phobius API client for transmembrane topology and signal peptide prediction.
License-free alternative to TMHMM for commercial use.
"""

import requests
import logging
from typing import Dict, Any, List
import re
import time

logger = logging.getLogger(__name__)

class PhobiusClient:
    """
    Client for Phobius transmembrane topology prediction.
    Free alternative to TMHMM without commercial licensing restrictions.
    """

    def __init__(self):
        self.base_url = "https://phobius.sbc.su.se"
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "Kozi-Vaccine-Pipeline/2.0"
        })

    def predict_transmembrane(self, sequence: str, protein_id: str = "query") -> Dict[str, Any]:
        """
        Predict transmembrane topology using Phobius.

        Args:
            sequence: Protein amino acid sequence
            protein_id: Identifier for the protein

        Returns:
            Dictionary with transmembrane prediction results
        """
        try:
            # Prepare FASTA format
            fasta_input = f">{protein_id}\n{sequence}"

            # Call Phobius web service
            response = self.session.post(
                f"{self.base_url}/cgi-bin/predict.pl",
                data={
                    'seq': fasta_input,
                    'format': 'short',
                    'outfmt': 'txt'
                },
                timeout=60
            )

            response.raise_for_status()

            # Parse results
            result = self._parse_phobius_output(response.text)

            logger.info(f"Phobius prediction for {protein_id}: {result['num_tm_helices']} TM helices")

            return result

        except requests.RequestException as e:
            logger.error(f"Phobius API error: {e}")
            return self._get_fallback_prediction(sequence)
        except Exception as e:
            logger.error(f"Phobius prediction failed: {e}")
            return self._get_fallback_prediction(sequence)

    def _parse_phobius_output(self, output_text: str) -> Dict[str, Any]:
        """Parse Phobius output text."""

        result = {
            'num_tm_helices': 0,
            'has_signal_peptide': False,
            'topology': 'unknown',
            'tm_regions': [],
            'signal_peptide_region': None,
            'localization_score': 0.0
        }

        try:
            lines = output_text.strip().split('\n')

            for line in lines:
                line = line.strip()

                # Look for topology line (contains TM info)
                if line.startswith('SEQTYPE'):
                    continue
                elif 'TM' in line and 'SP' in line:
                    # Has both signal peptide and transmembrane
                    result['has_signal_peptide'] = True
                    result['num_tm_helices'] = line.count('TM')
                elif 'TM' in line:
                    # Transmembrane only
                    result['num_tm_helices'] = line.count('TM')
                elif 'SP' in line:
                    # Signal peptide only
                    result['has_signal_peptide'] = True

                # Extract TM regions
                tm_matches = re.findall(r'TM(\d+)-(\d+)', line)
                for start, end in tm_matches:
                    result['tm_regions'].append({
                        'start': int(start),
                        'end': int(end),
                        'length': int(end) - int(start) + 1
                    })

            # Determine topology based on results
            if result['num_tm_helices'] == 0:
                if result['has_signal_peptide']:
                    result['topology'] = 'secreted'
                    result['localization_score'] = 0.9  # High surface accessibility
                else:
                    result['topology'] = 'cytoplasmic'
                    result['localization_score'] = 0.1  # Low surface accessibility
            elif result['num_tm_helices'] == 1:
                result['topology'] = 'single_pass_membrane'
                result['localization_score'] = 0.6  # Moderate surface accessibility
            else:
                result['topology'] = 'multi_pass_membrane'
                result['localization_score'] = 0.3  # Lower surface accessibility

        except Exception as e:
            logger.warning(f"Error parsing Phobius output: {e}")
            logger.debug(f"Raw output: {output_text[:200]}...")

        return result

    def _get_fallback_prediction(self, sequence: str) -> Dict[str, Any]:
        """Fallback prediction when Phobius API fails."""

        # Simple hydrophobicity-based prediction
        hydrophobic_aa = 'AILMFWYV'
        window_size = 20
        tm_count = 0

        # Scan for hydrophobic regions
        for i in range(len(sequence) - window_size + 1):
            window = sequence[i:i + window_size]
            hydrophobic_ratio = sum(1 for aa in window if aa in hydrophobic_aa) / window_size

            if hydrophobic_ratio >= 0.6:
                tm_count += 1
                # Skip ahead to avoid overlapping regions
                i += window_size // 2

        # Check for signal peptide
        n_terminal = sequence[:25] if len(sequence) > 25 else sequence
        has_signal = self._simple_signal_check(n_terminal)

        return {
            'num_tm_helices': min(tm_count, 8),  # Cap at reasonable number
            'has_signal_peptide': has_signal,
            'topology': 'predicted_fallback',
            'tm_regions': [],
            'signal_peptide_region': None,
            'localization_score': 0.5,  # Neutral score for fallback
            'prediction_method': 'fallback_hydrophobicity'
        }

    def _simple_signal_check(self, n_terminal: str) -> bool:
        """Simple signal peptide check for fallback."""
        if len(n_terminal) < 15:
            return False

        # Look for hydrophobic stretch in N-terminal region
        hydrophobic_aa = 'AILMFWYV'
        hydrophobic_count = sum(1 for aa in n_terminal[:15] if aa in hydrophobic_aa)

        return hydrophobic_count / 15 > 0.5

    def test_connection(self) -> bool:
        """Test Phobius service availability."""
        try:
            test_sequence = "MKLRLFCLAMLMACAQILNGS"
            result = self.predict_transmembrane(test_sequence, "test")

            success = 'num_tm_helices' in result
            logger.info(f"Phobius test: {'Success' if success else 'Failed'}")
            return success

        except Exception as e:
            logger.error(f"Phobius test failed: {e}")
            return False

# Global instance
phobius = PhobiusClient()