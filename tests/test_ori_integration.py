from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from plasmid_builder.assembler import assemble_plasmid


def test_ori_is_inserted():
    backbone = SeqRecord(Seq("AAAA"), id="backbone")
    backbone.features = []

    ori_info = {
        "ori_sequence": "TTTT",
        "method": "test"
    }

    antibiotic_markers = []
    mcs_seq = ""
    insert = SeqRecord(Seq("CCCC"), id="insert")

    plasmid = assemble_plasmid(
        backbone_record=backbone,
        ori_info=ori_info,
        antibiotic_markers=antibiotic_markers,
        mcs_sequence=mcs_seq,
        insert_record=insert
    )

    assert "TTTT" in str(plasmid.seq)
    assert any(f.type == "rep_origin" for f in plasmid.features)
