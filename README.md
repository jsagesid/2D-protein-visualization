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
