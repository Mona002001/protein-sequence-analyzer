import pytest
from analyzer import analyze_protein, validate_sequence
from pathlib import Path
from report_generator import generate_tsv, generate_html, generate_pdf


def test_sequence_validation():
    assert validate_sequence("ACDEFG") == "ACDEFG"


def test_invalid_sequence():
    with pytest.raises(ValueError):
        validate_sequence("ACDXYZ")


def test_protein_length():
    result = analyze_protein("ACDEFG")
    assert result["length"] == 6


def test_amino_acid_count():
    result = analyze_protein("ACDEFG")
    assert result["composition"]["A"]["count"] == 1


def test_amino_acid_percentages():
    result = analyze_protein("AAAA")
    assert result["composition"]["A"]["percentage"] == 100


def create_sample_report_data():
    """Create sample data for report generation tests."""
    results = analyze_protein("ACDEFGHIKLMNPQRSTVWY")

    return [{
        "protein_id": "Test_Protein",
        "results": results
    }]


def test_generate_tsv(tmp_path):
    records = create_sample_report_data()
    output_base = tmp_path / "test_report"

    output_path = generate_tsv(records, output_base)

    assert output_path.exists()
    assert output_path.suffix == ".tsv"
    assert "Test_Protein" in output_path.read_text()


def test_generate_html(tmp_path):
    records = create_sample_report_data()
    output_base = tmp_path / "test_report"

    output_path = generate_html(records, output_base)

    assert output_path.exists()
    assert output_path.suffix == ".html"

    content = output_path.read_text(encoding="utf-8")
    assert "Test_Protein" in content
    assert "Amino Acid Composition" in content


def test_generate_pdf(tmp_path):
    records = create_sample_report_data()
    output_base = tmp_path / "test_report"

    output_path = generate_pdf(records, output_base)

    assert output_path.exists()
    assert output_path.suffix == ".pdf"
    assert output_path.read_bytes().startswith(b"%PDF")