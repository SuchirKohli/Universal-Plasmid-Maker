# plasmid_builder/io.py
from typing import List, Tuple
import re

# The input can either be the name of a restriction enzyme (e.g., "EcoRI") or a recognition site (e.g., "GAATTC").

def parse_design(path: str) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str]]]:
    """
    Parse Design.txt.

    Returns:
        mcs_entries: List of (mcs_name, enzyme_or_site)
        antibiotic_markers: List of (marker_name, antibiotic_name)

    Accepts:
        MCS1, EcoRI
        MCS2, GAATTC
        kanR, Kanamycin
    """
    mcs_entries = []
    antibiotic_markers = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("**"):
                continue

            parts = [p.strip() for p in line.split(",", 1)]
            if len(parts) != 2:
                raise ValueError(f"Invalid format in Design.txt: {line}")

            left, right = parts

            # If right side looks like a DNA sequence → MCS with recognition site
            if re.fullmatch(r"[ACGTacgt]+", right):
                mcs_entries.append((left, right.upper()))
                continue

            # Heuristic: MCS lines usually contain these keywords
            if re.search(r"mcs|cloning|site|multiple", left, re.IGNORECASE):
                mcs_entries.append((left, right))
            else:
                antibiotic_markers.append((left, right))

    return mcs_entries, antibiotic_markers

def normalize_antibiotic_markers(antibiotic_entries):
    """
    antibiotic_entries: list of (sequence, antibiotic_name)
    """
    markers = []
    for seq, antibiotic in antibiotic_entries:
        seq = seq.strip().upper()
        if not all(base in "ACGT" for base in seq):
            raise ValueError(f"Invalid DNA sequence for antibiotic marker: {seq[:20]}...")
        
        markers.append({
            "sequence": seq,
            "antibiotic": antibiotic.strip()
        })
    return markers
