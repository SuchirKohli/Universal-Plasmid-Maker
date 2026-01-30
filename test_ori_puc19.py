from Bio import SeqIO
from plasmid_builder.ori_finder import predict_ori

# Path to pUC19 FASTA
fasta_path = (
    "/home/suchir/Documents/Suchir/Programming/CollegeProjects/"
    "BBL434 Lab/Assignment1/Universal-Plasmid-Maker/data/input/pUC19.fa"
)

# Load sequence
record = SeqIO.read(fasta_path, "fasta")
sequence = str(record.seq)

print(f"Sequence length: {len(sequence)} bp")

# Run ORI finder
ori = predict_ori(sequence)

print("\n=== ORI Prediction ===")
print(f"ORI start: {ori['ori_start']}")
print(f"ORI end:   {ori['ori_end']}")
print(f"ORI length: {ori['ori_end'] - ori['ori_start']} bp")
print("\nORI sequence (first 100 bp):")
print(ori['ori_sequence'][:100])
