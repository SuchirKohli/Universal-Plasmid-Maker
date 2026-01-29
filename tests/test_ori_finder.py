from plasmid_builder.ori_finder import predict_ori

def test_predict_ori_basic():
    seq = "A" * 5000 + "GCGCGCGCGC" * 100 + "A" * 5000
    ori = predict_ori(seq)

    assert "ori_sequence" in ori
    assert len(ori["ori_sequence"]) > 0
    assert ori["ori_start"] < ori["ori_end"]
