from plasmid_builder.io import parse_design
from plasmid_builder.io import normalize_antibiotic_markers
import pytest


def test_parse_design_mixed_entries(tmp_path):
    design = tmp_path / "Design.txt"
    design.write_text(
        """
        Multiple_Cloning_Site1, EcoRI
        Multiple_Cloning_Site2, GAATTC
        kanR, Kanamycin
        ampR, Ampicillin
        ** End of File
        """
    )

    mcs, antibiotics = parse_design(design)

    assert len(mcs) == 2
    assert mcs[0] == ("Multiple_Cloning_Site1", "EcoRI")
    assert mcs[1] == ("Multiple_Cloning_Site2", "GAATTC")

    assert len(antibiotics) == 2
    assert antibiotics[0] == ("kanR", "Kanamycin")
    assert antibiotics[1] == ("ampR", "Ampicillin")

def test_normalize_antibiotic_markers():
    entries = [
        ("atgcccggg", "Kanamycin")
    ]

    result = normalize_antibiotic_markers(entries)

    assert result[0]["sequence"] == "ATGCCCGGG"
    assert result[0]["antibiotic"] == "Kanamycin"


def test_invalid_antibiotic_sequence():
    entries = [
        ("ATGXYZ", "Kanamycin")
    ]

    with pytest.raises(ValueError):
        normalize_antibiotic_markers(entries)
