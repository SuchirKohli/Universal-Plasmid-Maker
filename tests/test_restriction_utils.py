from plasmid_builder.restriction_utils import (
    resolve_enzyme_or_site,
    resolve_mcs_entries
)
import pytest


def test_resolve_enzyme_name():
    result = resolve_enzyme_or_site("EcoRI")
    assert result["type"] == "enzyme"
    assert result["site"] == "GAATTC"


def test_resolve_recognition_site():
    result = resolve_enzyme_or_site("GAATTC")
    assert result["type"] == "site"
    assert result["site"] == "GAATTC"


def test_invalid_enzyme():
    with pytest.raises(ValueError):
        resolve_enzyme_or_site("NotARealEnzyme")


def test_resolve_mcs_entries():
    mcs = [
        ("Multiple_Cloning_Site1", "EcoRI"),
        ("Multiple_Cloning_Site2", "GAATTC")
    ]
    resolved = resolve_mcs_entries(mcs)

    assert resolved[0]["site"] == "GAATTC"
    assert resolved[1]["type"] == "site"
