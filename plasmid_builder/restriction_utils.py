# plasmid_builder/restriction_utils.py
import re
from Bio import Restriction


def is_recognition_sequence(token: str) -> bool:
    """Return True if token looks like a DNA recognition site."""
    return bool(re.fullmatch(r"[ACGTacgt]+", token))


def resolve_enzyme_or_site(token: str) -> dict:
    """
    Resolve a restriction enzyme name or a raw recognition site.

    Returns a dictionary with keys:
        type: 'enzyme' or 'site'
        enzyme: enzyme name or None
        site: recognition sequence (uppercase)

    Raises:
        ValueError if enzyme name is invalid
    """
    token = token.strip()

    # Case 1: User directly provided recognition site (GAATTC)
    if is_recognition_sequence(token):
        return {
            "type": "site",
            "enzyme": None,
            "site": token.upper()
        }

    # Case 2: User provided enzyme name (EcoRI)
    try:
        enzyme_class = getattr(Restriction, token)
    except AttributeError:
        raise ValueError(f"Unknown restriction enzyme: {token}")

    return {
        "type": "enzyme",
        "enzyme": token,
        "site": str(enzyme_class.site)
    }

def resolve_mcs_entries(mcs_entries):
    """
    Takes output of parse_design() MCS list and resolves enzymes/sites.
    """
    resolved = []
    for name, token in mcs_entries:
        info = resolve_enzyme_or_site(token)
        info["name"] = name
        resolved.append(info)
    return resolved
