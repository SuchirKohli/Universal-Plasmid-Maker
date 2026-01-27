from plasmid_builder.backbone import load_default_backbone, annotate_replication_features


def test_annotate_rsf1010_features():
    record, metadata = load_default_backbone("rsf1010")
    record = annotate_replication_features(record, metadata)
    assert len(record.features) == 3
