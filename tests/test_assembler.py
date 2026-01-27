from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from plasmid_builder.assembler import assemble_plasmid


def test_assemble_plasmid_basic():
    backbone = SeqRecord(Seq("AAAA"), id="backbone")
    backbone.features = []

    antibiotic_markers = [
        {"sequence": "TTTT", "antibiotic": "Kanamycin"}
    ]

    mcs_seq = "GGGG"
    insert = SeqRecord(Seq("CCCC"), id="insert")

    plasmid = assemble_plasmid(
        backbone,
        antibiotic_markers,
        mcs_seq,
        insert
    )

    assert str(plasmid.seq) == "AAAATTTTGGGGCCCC"
    assert plasmid.annotations["topology"] == "circular"
    assert len(plasmid.features) == 3
