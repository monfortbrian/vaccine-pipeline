import logging
from typing import Dict, Any, Optional

logger = logging.getLogger(__name__)

class PSORTbClient:
    def __init__(self):
        pass

    def predict_localization(self, sequence: str, organism_type: str = "gram_negative") -> Dict[str, Any]:
        # Try real PSORTb Docker first
        try:
            result = self._run_psortb_docker(sequence, organism_type)
            if result and result.get("localization") != "unknown":
                logger.info(f"PSORTb Docker: {result['localization']}")
                return result
        except Exception as e:
            logger.debug(f"PSORTb Docker not available: {e}")

        # Fallback to rule-based
        try:
            return self._predict_localization_rules(sequence, organism_type)
        except Exception as e:
            logger.error(f"PSORTb prediction failed: {e}")
            return {"localization": "unknown", "status": "error"}

    def _run_psortb_docker(self, sequence: str, organism_type: str) -> Optional[Dict[str, Any]]:
        """Run real PSORTb v3.0 via Docker container."""
        import subprocess
        import tempfile
        import os

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False, dir='.') as f:
            f.write(f">query\n{sequence}\n")
            fasta_path = f.name

        try:
            gram = "-n" if organism_type == "gram_negative" else "-p"
            abs_path = os.path.abspath(fasta_path)
            abs_dir = os.path.dirname(abs_path)
            filename = os.path.basename(abs_path)

            cmd = [
                "docker", "run", "--rm",
                "-v", f"{abs_dir}:/input",
                "brinkmanlab/psortb_commandline:1.0.2",
                gram, f"/input/{filename}"
            ]

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)

            if result.returncode != 0:
                logger.warning(f"PSORTb Docker error: {result.stderr[:200]}")
                return None

            return self._parse_psortb_output(result.stdout)
        finally:
            try:
                os.unlink(fasta_path)
            except Exception:
                pass

    def _parse_psortb_output(self, output: str) -> Optional[Dict[str, Any]]:
        """Parse PSORTb v3.0 text output."""
        if not output.strip():
            return None

        localization = "unknown"
        confidence = 0.0
        scores = {}

        for line in output.strip().split('\n'):
            line = line.strip()
            if 'final' in line.lower() and ('prediction' in line.lower() or 'localization' in line.lower()):
                loc_text = line.split(':')[-1].strip().lower() if ':' in line else line.split()[-1].lower()
                if 'cytoplasmic' in loc_text and 'membrane' not in loc_text:
                    localization = "cytoplasmic"
                elif 'cytoplasmic membrane' in loc_text or 'inner membrane' in loc_text:
                    localization = "inner_membrane"
                elif 'periplasmic' in loc_text:
                    localization = "periplasmic"
                elif 'outer membrane' in loc_text:
                    localization = "outer_membrane"
                elif 'extracellular' in loc_text:
                    localization = "extracellular"
                elif 'unknown' in loc_text:
                    localization = "unknown"
            else:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        score_val = float(parts[-1])
                        loc_name = ' '.join(parts[:-1]).lower()
                        scores[loc_name] = score_val
                        if score_val > confidence:
                            confidence = score_val
                    except ValueError:
                        pass

        surface_map = {
            "extracellular": 1.0, "outer_membrane": 0.9,
            "periplasmic": 0.6, "inner_membrane": 0.3,
            "cytoplasmic": 0.1, "unknown": 0.0,
        }

        return {
            "localization": localization,
            "confidence": min(confidence / 10.0, 1.0),
            "surface_score": surface_map.get(localization, 0.0),
            "scores": scores,
            "features": {"method": "PSORTb_3.0_Docker"},
        }

    def _predict_localization_rules(self, sequence: str, organism_type: str) -> Dict[str, Any]:
        """
        Rule-based localization prediction using signal peptides and transmembrane regions.
        """
        n_terminal = sequence[:30] if len(sequence) > 30 else sequence
        has_signal = self._check_signal_peptide(n_terminal)
        tm_regions = self._count_transmembrane_regions(sequence)
        has_lipoprotein_signal = self._check_lipoprotein_signal(sequence)

        if has_lipoprotein_signal:
            localization = "outer_membrane"
            confidence = 0.8
        elif tm_regions > 2:
            localization = "inner_membrane"
            confidence = 0.7
        elif has_signal and tm_regions == 0:
            if organism_type == "gram_positive":
                localization = "extracellular"
            else:
                localization = "periplasmic"
            confidence = 0.6
        elif tm_regions == 1:
            localization = "inner_membrane"
            confidence = 0.5
        else:
            localization = "cytoplasmic"
            confidence = 0.7

        if localization in ["outer_membrane", "extracellular"]:
            surface_score = 1.0
        elif localization == "periplasmic":
            surface_score = 0.6
        elif localization == "inner_membrane":
            surface_score = 0.3
        else:
            surface_score = 0.1

        return {
            "localization": localization,
            "confidence": confidence,
            "surface_score": surface_score,
            "scores": {
                "cytoplasmic": confidence if localization == "cytoplasmic" else 0.2,
                "inner_membrane": confidence if localization == "inner_membrane" else 0.1,
                "periplasmic": confidence if localization == "periplasmic" else 0.1,
                "outer_membrane": confidence if localization == "outer_membrane" else 0.1,
                "extracellular": confidence if localization == "extracellular" else 0.1
            },
            "features": {
                "signal_peptide": has_signal,
                "transmembrane_regions": tm_regions,
                "lipoprotein_signal": has_lipoprotein_signal
            }
        }

    def _check_signal_peptide(self, n_terminal: str) -> bool:
        if len(n_terminal) < 15:
            return False
        n_region = n_terminal[:5]
        h_region = n_terminal[5:15] if len(n_terminal) >= 15 else n_terminal[5:]
        positive_charges = sum(1 for aa in n_region if aa in 'KR')
        hydrophobic = sum(1 for aa in h_region if aa in 'AILMFWYV')
        hydrophobic_ratio = hydrophobic / len(h_region) if h_region else 0
        return positive_charges >= 1 and hydrophobic_ratio > 0.4

    def _count_transmembrane_regions(self, sequence: str) -> int:
        hydrophobic_aa = 'AILMFWYV'
        tm_count = 0
        window_size = 20
        hydrophobic_threshold = 0.65
        i = 0
        while i < len(sequence) - window_size + 1:
            window = sequence[i:i + window_size]
            hydrophobic_ratio = sum(1 for aa in window if aa in hydrophobic_aa) / window_size
            if hydrophobic_ratio >= hydrophobic_threshold:
                tm_count += 1
                i += window_size
            else:
                i += 1
        return tm_count

    def _check_lipoprotein_signal(self, sequence: str) -> bool:
        if len(sequence) < 20:
            return False
        n_terminal = sequence[:20]
        for i in range(1, len(n_terminal)):
            if n_terminal[i] == 'C':
                if i >= 2:
                    pattern = n_terminal[i-2:i+1]
                    if (pattern[0] in 'LVI' and
                        pattern[1] in 'ASTVI' and
                        pattern[2] == 'C'):
                        return True
        return False

# Global instance
psortb = PSORTbClient()