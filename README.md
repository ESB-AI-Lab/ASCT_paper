# Comparative Genomics of the ASCT Gene Family Across Mollusca

Code and analysis pipelines for the manuscript on Alanine, Serine, Cysteine Transporter (ASCT) gene family evolution across molluscan genomes.

## Overview

This repository contains reproducible workflows for identifying ASCT gene copies across 520 molluscan and outgroup genomes, reconstructing their phylogenetic relationships, and visualizing synteny and copy number variation. The focal taxa are six Mytilidae species (*Mytilus californianus*, *M. galloprovincialis*, *M. edulis*, *M. coruscus*, *Perna viridis*, and *Arcuatula senhousia*), with 15 additional outgroup species spanning Bivalvia, Gastropoda, Cephalopoda, Polyplacophora, and Annelida.

## Repository Structure

```
.
├── alignment_n_extraction/     # Gene discovery and region extraction
│   ├── sbatch_*.sh             # SLURM batch scripts for BLAST searches
│   ├── bash_*.sh               # Shell pipelines for BED/FASTA processing
│   ├── identify_gene_copies.py # Group BLAST hits into gene-copy loci
│   └── bed_files/              # Synteny figure pipeline
│       ├── config.yaml         # Central configuration for all synteny scripts
│       ├── extract_target_contigs.py
│       ├── 00_preprocess.py
│       ├── 01_synteny_plot.py
│       └── 02_composite-figure.py
├── phylogenetics/              # Alignment, tree inference, and comparison plots
│   ├── bash_prank.sh           # Multiple sequence alignment (PRANK)
│   ├── bash_mrbayes.sh         # Bayesian phylogenetic inference (MrBayes)
│   ├── biophylo.v2.py          # Consensus tree visualization
│   ├── gene_arrangement_plot.py
│   ├── gene_comparison.py      # Pairwise identity heatmap and ribbon plot
│   └── fastas/                 # Gene sequences and summary data
└── igv/                        # IGV.js genome browser notebooks
```

## Analysis Pipeline

### 1. Genome Survey and Gene Discovery

BLAST-based identification of ASCT homologs across molluscan genomes using the *Fasciola hepatica* ASCT protein (UniProt C6EUD4) as query.

```bash
# Download genomes from NCBI (SLURM array job, 520 accessions)
sbatch sbatch_datasets.sh

# Build BLAST databases
sbatch sbatch_blast_buildgenome.sh

# Run tBLASTn (protein query vs. nucleotide databases)
sbatch sbatch_tblastn.sh

# Run BLASTn (nucleotide query vs. nucleotide databases)
sbatch sbatch_blast.sh
```

### 2. Hit Processing and Copy Identification

```bash
# Convert BLAST outfmt6 to BED format with strand-aware merging
bash bash_outfmt6_to_bed.sh <reference_prefix>

# Generate per-species ASCT BED files with copy numbering
bash bash_identify_gene_copies.sh

# Extract exon sequences from genomes
bash bash_extract_fastas.sh

# Concatenate exons into gene-level FASTAs
cd bed_files/
bash bash_extract_genes.sh
```

### 3. Transposable Element Annotation

```bash
# Extract flanking regions and run RepeatMasker
bash bash_extract_mytilus_regions.sh

# Generate TE landscape CSV
bash bash_generate_ASCT_TE_annotation_csv.sh
```

### 4. Phylogenetic Analysis

```bash
cd phylogenetics/

# Multiple sequence alignment
bash bash_prank.sh all_ASCT_genes.fasta

# Bayesian inference (GTR+I+G, 1M generations, 4 chains x 2 runs)
bash bash_mrbayes.sh

# Visualize consensus tree
python biophylo.v2.py

# All-vs-all pairwise identity (MAFFT) + heatmap + ribbon plot
python gene_comparison.py

# Gene arrangement figure
python gene_arrangement_plot.py
```

### 5. Synteny Figure

Config-driven pipeline producing the main manuscript figure. Edit `bed_files/config.yaml` to adjust species, alignment parameters, and plot settings.

```bash
cd alignment_n_extraction/bed_files/

# Extract target contigs/regions from genomes
python extract_target_contigs.py --config config.yaml

# Parse species tree + run pairwise alignments (minimap2)
python 00_preprocess.py --config config.yaml

# Generate synteny diagram (Panel A) and presence/absence heatmap (Panel B)
python 01_synteny_plot.py --config config.yaml

# Assemble composite figure + LaTeX template
python 02_composite-figure.py --config config.yaml
```

## Dependencies

### HPC / Command-line Tools

| Tool | Version | Purpose |
|------|---------|---------|
| [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/) | - | Genome download |
| [BLAST+](https://blast.ncbi.nlm.nih.gov/) | 2.2.31+ | Sequence homology search |
| [samtools](http://www.htslib.org/) | - | FASTA indexing and region extraction |
| [bedtools](https://bedtools.readthedocs.io/) | - | BED manipulation and sequence extraction |
| [PRANK](http://wasabiapp.org/software/prank/) | - | Phylogeny-aware multiple alignment |
| [MrBayes](https://nbisweden.github.io/MrBayes/) | - | Bayesian phylogenetic inference |
| [RepeatMasker](http://www.repeatmasker.org/) | - | Transposable element annotation |
| [minimap2](https://github.com/lh3/minimap2) | 2.30 | Whole-genome pairwise alignment |
| [MAFFT](https://mafft.cbrc.jp/alignment/software/) | - | Pairwise sequence alignment |

### Python (3.9+)

```
biopython
matplotlib
seaborn
numpy
pandas
scipy
pyyaml
pygenomeviz
```

## Species Codes

Four-letter abbreviations used across all scripts, BED files, and figures:

| Code | Species | Common Name |
|------|---------|-------------|
| Mcal | *Mytilus californianus* | California mussel |
| Mgal | *Mytilus galloprovincialis* | Mediterranean mussel |
| Mcor | *Mytilus coruscus* | Korean mussel |
| Medi | *Mytilus edulis* | Blue mussel |
| Pvir | *Perna viridis* | Asian green mussel |
| Asen | *Arcuatula senhousia* | Asian date mussel |

Outgroup species: Fhep (*Fasciola hepatica*), Svel (*Solemya velum*), Pvul (*Patella vulgata*), Hruf (*Haliotis rufescens*), Pcan (*Pomacea canaliculata*), Osin (*Octopus sinensis*), Npom (*Nautilus pompilius*), Lsin (*Liolophura sinensis*), Ofus (*Owenia fusiformis*), Mmre (*Meretrix meretrix*), Mmar (*Margaritifera margaritifera*), Upic (*Unio pictorum*), Mmrc (*Mercenaria mercenaria*), Rphi (*Ruditapes philippinarum*), Pmax (*Pecten maximus*), Mgig (*Magallana gigas*), Airr (*Argopecten irradians*).

## Data Availability

Genome assemblies were obtained from NCBI GenBank/RefSeq. Accession numbers for all 520 genomes surveyed are listed in `alignment_n_extraction/dbs_dl.txt`. The protein query sequence corresponds to *Fasciola hepatica* ASCT (UniProt [C6EUD4](https://www.uniprot.org/uniprot/C6EUD4)).

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) for details.

Copyright (c) 2026 ESB AI Lab
