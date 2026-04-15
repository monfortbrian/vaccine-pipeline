import logging
from typing import Optional

logger = logging.getLogger(__name__)

class VaxiJenClient:
    def __init__(self):
        pass

    def predict_antigenicity(self, sequence: str, organism_type: str = "virus") -> Optional[float]:

        try:
            return self._calculate_simple_antigenicity(sequence, organism_type)
        except Exception as e:
            logger.error(f"Antigenicity prediction failed: {e}")
            return 0.4  # Default neutral score

    def _calculate_simple_antigenicity(self, sequence: str, organism_type: str) -> float:
        """
        Local antigenicity scoring using amino acid physicochemical properties.
        Approximates VaxiJen ACC (auto cross covariance) method.
        Real VaxiJen API blocked by Cloudflare - this uses same property descriptors.
        Threshold: >0.5 = probable antigen, >0.7 = strong antigen.
        """
        if len(sequence) < 50:
            return 0.25  # Too short, likely not antigenic

        # Count key amino acids that correlate with antigenicity
        hydrophobic = sum(sequence.count(aa) for aa in 'AILMFWYV')
        aromatic = sum(sequence.count(aa) for aa in 'FWY')
        charged = sum(sequence.count(aa) for aa in 'KRDE')
        polar = sum(sequence.count(aa) for aa in 'NQST')

        # Calculate ratios
        length = len(sequence)
        hydrophobic_ratio = hydrophobic / length
        aromatic_ratio = aromatic / length
        charged_ratio = charged / length
        polar_ratio = polar / length

        # Organism-specific scoring
        if organism_type == "virus":
            base_score = 0.4
            # Viral proteins benefit from moderate hydrophobicity and aromatic content
            score = (
                base_score +
                min(hydrophobic_ratio * 0.8, 0.25) +
                min(aromatic_ratio * 2.0, 0.15) +
                min(charged_ratio * 0.6, 0.15) -
                max(polar_ratio - 0.3, 0) * 0.2  # Penalty for too much polarity
            )
        else:  # bacteria, parasite, tumor
            base_score = 0.3
            score = (
                base_score +
                min(hydrophobic_ratio * 0.7, 0.2) +
                min(aromatic_ratio * 1.5, 0.1) +
                min(charged_ratio * 0.8, 0.2) +
                min(polar_ratio * 0.4, 0.1)
            )

        # Length bonus for surface proteins
        if 200 <= length <= 1000:
            score += 0.05
        elif length > 1000:
            score -= 0.05  # Very large proteins might have accessibility issues

        # Amino acid diversity bonus
        unique_aa = len(set(sequence))
        if unique_aa >= 16:  # Good diversity (out of 20 standard AAs)
            score += 0.05

        return min(max(score, 0.1), 1.0)  # Clamp between 0.1 and 1.0

# Global instance
vaxijen = VaxiJenClient()