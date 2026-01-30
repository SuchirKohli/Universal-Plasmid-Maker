from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from plasmid_builder.assembler import assemble_plasmid


def test_assemble_plasmid_basic():
    backbone = SeqRecord(Seq("AAAA"), id="backbone")
    backbone.features = []

    ori_info = {
        "ori_sequence": "TTTT",
        "method": "test"
    }

    antibiotic_markers = [
        {"sequence": "GGGG", "antibiotic": "Kanamycin"}
    ]

    mcs_seq = "CCCC"

    plasmid = assemble_plasmid(
        backbone_record=backbone,
        ori_info=ori_info,
        antibiotic_markers=antibiotic_markers,
        mcs_sequence=mcs_seq,
        insert_record=None
    )

    seq = str(plasmid.seq)
    assert "AAAA" in seq   # backbone
    assert "TTTT" in seq   # ORI
    assert "GGGG" in seq   # marker
    assert "CCCC" in seq   # MCS
