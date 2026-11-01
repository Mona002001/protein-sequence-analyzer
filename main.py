
import argparse
from pathlib import Path

from Bio import SeqIO
from analyzer import analyze_protein
from report_generator import generate_tsv, generate_html, generate_pdf


AVAILABLE_OPERATIONS = [
    "length",
    "composition",
    "molecular_weight",
    "isoelectric_point",
    "gravy",
    "instability_index",
    "aromaticity"
]


def read_fasta(input_file):
    """Read protein sequences from a FASTA file."""

    input_path = Path(input_file)

    if not input_path.exists():
        raise FileNotFoundError(
            f"Input file not found: {input_file}"
        )

    records = list(SeqIO.parse(input_path, "fasta"))

    if not records:
        raise ValueError("No sequences found in FASTA file.")

    return records


def display_results(record_id, results, operations):
    """Display selected analysis results."""

    print("\n" + "=" * 45)
    print(f"Protein ID: {record_id}")
    print("=" * 45)

    for operation in operations:

        if operation == "length":
            print(f"Sequence Length: {results['length']}")

        elif operation == "molecular_weight":
            print(
                f"Molecular Weight: "
                f"{results['molecular_weight']} Da"
            )

        elif operation == "isoelectric_point":
            print(
                f"Isoelectric Point: "
                f"{results['isoelectric_point']}"
            )

        elif operation == "gravy":
            print(f"GRAVY: {results['gravy']}")

        elif operation == "instability_index":
            print(
                f"Instability Index: "
                f"{results['instability_index']}"
            )

        elif operation == "aromaticity":
            print(f"Aromaticity: {results['aromaticity']}")

        elif operation == "composition":
            print("\nAmino Acid Composition:")

            for amino_acid, values in results["composition"].items():
                print(
                    f"{amino_acid}: "
                    f"{values['count']} "
                    f"({values['percentage']}%)"
                )


def main():

    parser = argparse.ArgumentParser(
        description="Protein Sequence Analyzer"
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input FASTA file"
    )

    parser.add_argument(
        "-o", "--output",
        default="output/protein_report",
        help="Output file name without extension"
    )

    parser.add_argument(
        "--operations",
        nargs="+",
        choices=AVAILABLE_OPERATIONS + ["all"],
        required=True,
        help="Select analysis operations"
    )

    args = parser.parse_args()

    if "all" in args.operations:
        operations = AVAILABLE_OPERATIONS
    else:
        operations = args.operations

    analyzed_records = []

    try:
        records = read_fasta(args.input)

        for record in records:

            sequence = str(record.seq)

            results = analyze_protein(sequence)

            analyzed_records.append({
                "protein_id": record.id,
                "results": results
            })

            display_results(
                record.id,
                results,
                operations
            )

        # Generate TSV report
        
        generate_tsv(analyzed_records, args.output)

        generate_html(analyzed_records, args.output)

        generate_pdf(analyzed_records, args.output)
        
        print("\nAnalysis completed successfully!")

        print(f"Output file base name: {args.output}")

    except (FileNotFoundError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()

