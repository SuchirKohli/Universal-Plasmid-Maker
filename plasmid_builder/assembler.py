# plasmid_builder/assembler.py

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


def assemble_plasmid(
    backbone_record,
    ori_info,
    antibiotic_markers,
    mcs_sequence,
    insert_record=None
):
    """
    Assemble final plasmid sequence.

    Parameters
    ----------
    backbone_record : SeqRecord
        RSF1010 backbone with replication genes annotated
    antibiotic_markers : list of dicts
        Output of normalize_antibiotic_markers()
    mcs_sequence : str
        Concatenated MCS sequence
    insert_record : SeqRecord or None
        Input insert DNA

    Returns
    -------
    SeqRecord
        Assembled plasmid
    """

    seq_parts = []
    features = []
    cursor = 0

    # Backbone
    seq_parts.append(backbone_record.seq)
    for feat in backbone_record.features:
        features.append(feat)
    cursor += len(backbone_record.seq)

    # ORI
    ori_seq = Seq(ori_info["ori_sequence"])
    ori_len = len(ori_seq)

    features.append(
        SeqFeature(
            FeatureLocation(cursor, cursor + ori_len),
            type="rep_origin",
            qualifiers={"note": "Predicted Origin of Replication",
                        "method": ori_info.get("method", "ORI finder")}
        )
    )
    seq_parts.append(ori_seq)
    cursor += ori_len

    # Antibiotic markers
    for marker in antibiotic_markers:
        seq = Seq(marker["sequence"])
        length = len(seq)

        features.append(
            SeqFeature(
                FeatureLocation(cursor, cursor + length),
                type="CDS",
                qualifiers={
                    "note": "antibiotic resistance marker",
                    "antibiotic": marker["antibiotic"]
                }
            )
        )

        seq_parts.append(seq)
        cursor += length

    # MCS
    mcs_seq = Seq(mcs_sequence)
    mcs_len = len(mcs_seq)

    features.append(
        SeqFeature(
            FeatureLocation(cursor, cursor + mcs_len),
            type="misc_feature",
            qualifiers={"note": "Multiple Cloning Site"}
        )
    )

    seq_parts.append(mcs_seq)
    cursor += mcs_len

    # Optional insert (not used in this assignment)
    if insert_record is not None:
        seq_parts.append(insert_record.seq)
        insert_len = len(insert_record.seq)

        features.append(
            SeqFeature(
                FeatureLocation(cursor, cursor + insert_len),
                type="misc_feature",
                qualifiers={"note": "Inserted sequence"}
            )
        )
        cursor += insert_len

    # Final plasmid
    final_seq = Seq("").join(seq_parts)

    plasmid = SeqRecord(
        final_seq,
        id="Assembled_Plasmid",
        description="Universal plasmid constructed in silico"
    )

    plasmid.features = features
    plasmid.annotations["topology"] = "circular"
    plasmid.annotations["molecule_type"] = "DNA"

    return plasmid

def build_mcs_sequence(resolved_mcs):
    """
    Concatenate recognition sites to form an MCS.
    """
    return "".join(entry["site"] for entry in resolved_mcs)
