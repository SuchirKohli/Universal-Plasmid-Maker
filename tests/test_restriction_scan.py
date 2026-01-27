from Bio.Seq import Seq
from plasmid_builder.restriction_utils import (
    scan_sequence_for_restriction_sites,
    detect_restriction_conflicts
)


def test_scan_known_enzyme():
    seq = Seq("AAAGAATTCTTT")
    resolved = [
        {
            "name": "Multiple_Cloning_Site1",
            "type": "enzyme",
            "enzyme": "EcoRI",
            "site": "GAATTC"
        }
    ]

    results = scan_sequence_for_restriction_sites(seq, resolved)
    assert results["Multiple_Cloning_Site1"] == [5]


def test_scan_raw_site():
    seq = Seq("CCCGGATCCAAA")
    resolved = [
        {
            "name": "Multiple_Cloning_Site2",
            "type": "site",
            "enzyme": None,
            "site": "GGATCC"
        }
    ]

    results = scan_sequence_for_restriction_sites(seq, resolved)
    assert results["Multiple_Cloning_Site2"] == [3]


def test_detect_conflicts():
    seq = Seq("AAAGAATTCTTT")
    resolved = [
        {
            "name": "MCS1",
            "type": "enzyme",
            "enzyme": "EcoRI",
            "site": "GAATTC"
        }
    ]

    conflicts = detect_restriction_conflicts(seq, resolved)
    assert "MCS1" in conflicts
