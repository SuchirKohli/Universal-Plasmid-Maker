# plasmid_builder/assembler.py

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


def assemble_plasmid(
    backbone_record,
    antibiotic_markers,
    mcs_sequence,
    insert_record
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
    insert_record : SeqRecord
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

    # Insert
    insert_len = len(insert_record.seq)

    features.append(
        SeqFeature(
            FeatureLocation(cursor, cursor + insert_len),
            type="insert",
            qualifiers={"note": "Inserted DNA"}
        )
    )

    seq_parts.append(insert_record.seq)
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

    return plasmid

def build_mcs_sequence(resolved_mcs):
    """
    Concatenate recognition sites to form an MCS.
    """
    return "".join(entry["site"] for entry in resolved_mcs)
