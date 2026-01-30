# plasmid_builder/restriction_utils.py
import re
from Bio import Restriction
from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch
import warnings

def is_recognition_sequence(token: str) -> bool:
    """Return True if token looks like a DNA recognition site."""
    return bool(re.fullmatch(r"[ACGTacgt]+", token))


import warnings
from Bio.Restriction import Restriction, AllEnzymes


def resolve_enzyme_or_site(token):
    """
    Resolve a restriction enzyme name or raw restriction site.
    Unknown enzymes raise a warning and return None.
    """

    token = token.strip()

    # Try raw DNA site
    if all(c in "ATGC" for c in token.upper()):
        return {
            "type": "site",
            "enzyme": None,
            "site": token.upper()
        }

    # Try enzyme name (case-insensitive)
    token_upper = token.upper()
    for enzyme in AllEnzymes:
        if enzyme.__name__.upper() == token_upper:
            return {
                "type": "enzyme",
                "enzyme": enzyme.__name__,
                "site": str(enzyme.site)
            }

    # If nothing worked → warn and skip
    warnings.warn(
        f"Unknown restriction enzyme '{token}'. This MCS entry will be skipped.",
        UserWarning
    )
    return None


def resolve_mcs_entries(mcs_entries):
    """
    Takes output of parse_design() MCS list and resolves enzymes/sites.
    """
    resolved = []
    for name, token in mcs_entries:
        info = resolve_enzyme_or_site(token)
        if info is not None:
            info["name"] = name
            resolved.append(info)
    return resolved

def scan_sequence_for_restriction_sites(seq: Seq, resolved_mcs):
    """
    Scan a DNA sequence for restriction enzyme cut sites.

    Parameters:
        seq : Bio.Seq.Seq
        resolved_mcs : list of resolved MCS dicts

    Returns:
        dict:
            key   -> MCS name
            value -> list of cut positions (0-based)
    """
    results = {}

    for entry in resolved_mcs:
        name = entry["name"]
        results[name] = []

        # Case 1: enzyme known → use Biopython
        if entry["type"] == "enzyme":
            enzyme_name = entry["enzyme"]
            enzyme = getattr(__import__("Bio.Restriction").Restriction, enzyme_name)
            rb = RestrictionBatch([enzyme])
            cuts = rb.search(seq)
            results[name] = cuts.get(enzyme, [])

        # Case 2: raw recognition site → string search
        else:
            site = entry["site"]
            site_len = len(site)
            for i in range(len(seq) - site_len + 1):
                if str(seq[i:i+site_len]).upper() == site:
                    results[name].append(i)

    return results

def detect_restriction_conflicts(seq: Seq, resolved_mcs):
    """
    Detect whether any restriction sites cut within a given sequence.

    Returns:
        dict:
            key   -> MCS name
            value -> number of cut sites found
    """
    scan_results = scan_sequence_for_restriction_sites(seq, resolved_mcs)
    conflicts = {}

    for name, positions in scan_results.items():
        if positions:
            conflicts[name] = positions

    return conflicts
