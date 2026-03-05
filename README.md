**2D Protein Visualization**

R tools for visualizing substitutions in focal species relative to annotated protein features and domains.

**Overview**

This repository contains an R script for visualizing amino-acid substitutions from multiple sequence alignments (MSAs) mapped onto known protein domains and functional features.
The workflow was designed for comparative genomics and evolutionary analyses, allowing users to examine how substitutions in focal species relate to known functional regions of a protein.

Key features of the script include:

- Mapping substitutions from aligned sequences onto protein feature coordinates<br>
- Correcting positional shifts caused by insertions and deletions (indels)<br>
- Importing protein feature annotations directly from UniProt<br>
- Generating 2D visualizations of proteins with annotated substitutions<br>
- Calculating summary statistics for amino acid diversity

This approach helps link sequence evolution with functional protein structure, which can be particularly useful when investigating candidate adaptive mutations.

**Input Files**

The script requires two main inputs.

1. Multiple Sequence Alignment

A trimmed protein alignment in FASTA format.

Example:

>Homo_sapiens
MTEYKLVVVGAGGVGKSALTIQLIQ
>Phoca_vitulina
MTEYKLVVVGAGGVGKSALTIQLIQ

Requirements:

- Protein alignment

- Trimmed to remove poorly aligned regions

- Sequence names must match those in the taxonomy file

2. Taxonomy File

A CSV file containing metadata for each sequence.

Example:

Species	seq.no	tax1	tax2	tax3<br>
human	1	Primates	Hominidae	Hominidae<br>
mouselemur	2	Primates	Cheirogaleidae	Cheirogaleidae<br>
mouse	3	Rodentia	Muridae	Muridae


Requirements:

- Must contain a column with sequence names

- Names must exactly match those used in the MSA

Output

The script generates a 2D visualization of the protein including:

- Annotated functional domains

- Highlighted substitutions in focal species

- Adjusted coordinates that account for alignment gaps
