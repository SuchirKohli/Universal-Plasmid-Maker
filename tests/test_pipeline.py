from plasmid_builder.pipeline import run_pipeline
from Bio import SeqIO
import os


def test_full_pipeline(tmp_path):
    input_fa = tmp_path / "input.fa"
    design = tmp_path / "Design.txt"
    output_fa = tmp_path / "Output.Fa"
    output_gb = tmp_path / "Output.gb"

    input_fa.write_text(">test\nATGCATGC")
    design.write_text("""
Multiple_Cloning_Site1, EcoRI
ATGCCC, Kanamycin
""")

    run_pipeline(
        input_fasta=str(input_fa),
        design_file=str(design),
        output_fasta=str(output_fa),
        output_genbank=str(output_gb)
    )

    assert os.path.exists(output_fa)
    record = SeqIO.read(output_fa, "fasta")
    assert len(record.seq) > 0
