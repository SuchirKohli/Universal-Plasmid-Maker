# Universal Plasmid Maker

A modular bioinformatics pipeline to construct a **universal plasmid** *in silico* from:
- a DNA insert sequence, and  
- a user-defined plasmid design specification.


---

## Overview

The pipeline takes two inputs:

1. **`Input.Fa`**  
   A FASTA file containing the DNA sequence of the organism in which the plasmid needs to be inserted. 

2. **`Design.txt`**  
   A text file specifying:
   - restriction enzymes or recognition sequences for constructing a Multiple Cloning Site (MCS),
   - antibiotic resistance marker(s), provided directly as DNA sequences along with the corresponding antibiotic name.

The pipeline produces a fully assembled plasmid sequence containing:
- origin of replication of the unknown organism,
- annotated replication genes,
- annotated antibiotic resistance marker(s),
- annotated MCS and insert regions,
- circular topology.

---

## Biological Design Choices

### Backbone Selection
- **RSF1010 (NCBI GenBank accession: M28829.1)** is used as the default plasmid backbone.
- RSF1010 is a **broad-host-range plasmid**, making it suitable for a “universal” plasmid design that is not limited to *Escherichia coli*–specific replication systems.

The complete nucleotide sequence of RSF1010 is sourced directly from **NCBI GenBank**, ensuring biological accuracy and traceability.

---

### Replication Genes
Replication machinery is annotated based on **curated GenBank CDS annotations** for RSF1010:
- `repA`
- `repB`
- `repC`

These genes constitute the minimal replication system required for autonomous plasmid replication across a broad range of bacterial hosts.

The origin of replication (`oriV`) is **not explicitly annotated as a discrete feature in the GenBank record** for RSF1010 and is therefore **not modeled as a separate annotation**. Replication is represented via the required replication genes, consistent with available database annotations. It can be observed from the research paper [4] that there also exists other essential genes (plasmid mobilization, oriV and oriT).

---

### Antibiotic Resistance Markers
- Antibiotic resistance markers are provided **directly as DNA sequences** in `Design.txt`.
- Each marker is paired with the name of the antibiotic it confers resistance to.

---

## Input Formats

### `Input.Fa`
```fasta
>insert_sequence
ATGCATGCATGC...
```
### `Design.txt`
```txt
Multiple_Cloning_Site1, EcoRI
Multiple_Cloning_Site2, GAATTC
ATGCCCAGTACG..., Kanamycin
```
## Environment Setup
```bash
conda env create -f environment.yml
conda activate plasmid-maker
```
---

## References
1. NCBI GenBank,
RSF1010 complete plasmid sequence,
Accession: M28829.1,
National Center for Biotechnology Information (NCBI)

2. Scholz, P., et al.,
Broad-host-range plasmids: replication and function,
FEMS Microbiology Reviews

3. Cock, P. J. A., et al.,
Biopython: freely available Python tools for computational molecular biology and bioinformatics,
Bioinformatics

4. Scholz, et al.,  Compiete nucleotide sequence and gene organization of the broad-host-range plasmid RSFlOlO (1988)