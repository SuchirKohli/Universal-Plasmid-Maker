from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from pathlib import Path
import json


BACKBONE_DIR = Path(__file__).parent / "backbones"


def load_default_backbone(name: str = "rsf1010"):
    """
    Load plasmid backbone sequence and metadata.

    Parameters
    ----------
    name : str
        Backbone name (folder name under backbones/)

    Returns
    -------
    record : SeqRecord
        Backbone sequence
    metadata : dict
        Backbone feature metadata
    """
    backbone_path = BACKBONE_DIR / name / f"{name}.fasta"
    metadata_path = BACKBONE_DIR / name / f"{name}.json"

    if not backbone_path.exists():
        raise FileNotFoundError(f"Backbone FASTA not found: {backbone_path}")

    if not metadata_path.exists():
        raise FileNotFoundError(f"Backbone metadata not found: {metadata_path}")

    record = SeqIO.read(backbone_path, "fasta")

    with open(metadata_path) as f:
        metadata = json.load(f)

    record.id = metadata.get("name", record.id)
    record.description = metadata.get("description", "")

    return record, metadata


def annotate_replication_features(record, metadata):
    """
    Annotate replication genes onto the backbone using metadata.
    """
    for gene in metadata.get("replication_genes", []):
        start = gene["start"]
        end = gene["end"]
        strand = gene.get("strand", 1)

        qualifiers = {
            "gene": gene["gene"],
            "product": gene.get("product", "")
        }

        if "protein_id" in gene:
            qualifiers["protein_id"] = gene["protein_id"]

        feature = SeqFeature(
            FeatureLocation(start, end, strand=strand),
            type="CDS",
            qualifiers=qualifiers
        )

        record.features.append(feature)

    return record
