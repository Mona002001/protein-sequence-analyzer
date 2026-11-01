
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Standard 20 amino acids
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def validate_sequence(sequence):
    """Validate and clean a protein sequence."""

    sequence = "".join(sequence.split()).upper()

    if not sequence:
        raise ValueError("Protein sequence is empty.")

    invalid = set(sequence) - set(AMINO_ACIDS)

    if invalid:
        raise ValueError(
            f"Invalid amino acids found: {', '.join(sorted(invalid))}"
        )

    return sequence


def analyze_protein(sequence):
    """Calculate physicochemical properties of a protein."""

    sequence = validate_sequence(sequence)

    protein = ProteinAnalysis(sequence)

    # Sequence length
    length = len(sequence)

    # Amino acid composition
    counts = protein.count_amino_acids()

    composition = {}

    for amino_acid in AMINO_ACIDS:
        count = counts.get(amino_acid, 0)

        percentage = (count / length) * 100

        composition[amino_acid] = {
            "count": count,
            "percentage": round(percentage, 2)
        }

    # Physicochemical properties
    results = {
        "sequence": sequence,
        "length": length,
        "molecular_weight": round(
            protein.molecular_weight(), 2
        ),
        "isoelectric_point": round(
            protein.isoelectric_point(), 2
        ),
        "gravy": round(
            protein.gravy(), 4
        ),
        "instability_index": round(
            protein.instability_index(), 2
        ),
        "aromaticity": round(
            protein.aromaticity(), 4
        ),
        "composition": composition
    }

    return results

