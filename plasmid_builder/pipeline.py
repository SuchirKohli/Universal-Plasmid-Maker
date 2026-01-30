from Bio import SeqIO
import os

from plasmid_builder.io import (
    parse_design,
    normalize_antibiotic_markers
)
from plasmid_builder.restriction_utils import resolve_mcs_entries
from plasmid_builder.backbone import load_default_backbone, annotate_replication_features
from plasmid_builder.assembler import assemble_plasmid, build_mcs_sequence
from plasmid_builder.ori_finder import predict_ori


def run_pipeline(
    input_fasta,
    design_file,
    output_fasta,
    output_genbank
):
    
    # Load host organism sequence
    host_record = SeqIO.read(input_fasta, "fasta")
    host_seq = str(host_record.seq)

    # Predict ORI from host
    ori_info = predict_ori(host_seq)


    # Parse Design.txt
    mcs_entries, antibiotic_entries = parse_design(design_file)

    # Normalize antibiotic markers
    antibiotic_markers = normalize_antibiotic_markers(antibiotic_entries)

    # Resolve MCS enzymes/sites
    resolved_mcs = resolve_mcs_entries(mcs_entries)

    # Build MCS sequence
    mcs_sequence = build_mcs_sequence(resolved_mcs)

    # Load and annotate backbone
    backbone_record, metadata = load_default_backbone("rsf1010")
    backbone_record = annotate_replication_features(backbone_record, metadata)

    # Assemble plasmid
    plasmid = assemble_plasmid(
        backbone_record,
        ori_info,
        antibiotic_markers,
        mcs_sequence,
        insert_record=None
    )

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_fasta), exist_ok=True)
    os.makedirs(os.path.dirname(output_genbank), exist_ok=True)

    # Write outputs
    SeqIO.write(plasmid, output_fasta, "fasta")
    SeqIO.write(plasmid, output_genbank, "genbank")

    return plasmid
