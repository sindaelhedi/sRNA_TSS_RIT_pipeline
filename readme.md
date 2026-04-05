# sRNA Discovery Pipeline

A comprehensive, generalizable bioinformatic pipeline for systematic discovery and validation of small regulatory RNAs (sRNAs) in bacterial genomes.

## Overview

This pipeline implements a sequential filtering strategy for bacterial sRNA discovery by integrating three independent biological constraints:

1. sRNA detection via sRNA-Detect
2. Transcription Start Site (TSS) Detection via dRNA-seq
3. Rho-Independent Terminator (RIT) Prediction via covariance models (RNIE)

 What This Pipeline Does

 Input
- **Genome sequence** (FASTA)
- **TSS predictions** (BED, from TSSAR)
- **sRNA candidates** (GTF, from sRNA-Detect)
- **BLAST reference database** (nucleotide)

 Output
- **Summary table** (TSV): sRNA candidates × score combinations with confidence metrics
- **Filtered BED files** (per threshold): sRNA, TSS, intersections
- **RNIE predictions** (GFF): Intrinsic terminators
- **BLAST results** (tabular): Homology validation
- **FASTA sequences**: For downstream analysis


 Installation

 Prerequisites

- **bedtools** ≥ 2.30
- **BLAST+** ≥ 2.14.1
- **SAMtools** ≥ 1.15
- **Perl** ≥ 5.30
- **RNIE** with covariance models

### Quick Setup

bash
# Clone repository
git clone https://github.com/sindalehedi/sRNA-discovery-pipeline.git
cd sRNA-discovery-pipeline


# Verify installation
bedtools --version
blastn -version
samtools --version
perl -v
```

### Install RNIE (if not already installed)

```bash
# Clone RNIE
git clone https://github.com/rnie/rnie.git /opt/RNIE
cd /opt/RNIE


## Usage

### Basic Command

```bash
bash pipeline.sh <BASE_DIR> <OUTPUT_ROOT> 
```

## Input Directory Structure

The pipeline expects the following structure:

```
BASE_DIR/
├── organism1/
│   ├── genome.fa              (or .fasta, .fna)
│   ├── tss.bed                (BED6 format, from TSSAR)
│   ├── srna.gtf               (GTF format, from sRNA-Detect)
│   └── blast_db/
│       ├── organism.nin       (BLAST index files)
│       ├── organism.nhr
│       └── organism.nsq
├── organism2/
│   ├── genome.fa
│   ├── tss.bed
│   ├── srna.gtf
│   └── blast_db/
└── ...


## Output Directory Structure

```
OUTPUT_ROOT/
├── organism1/
│   ├── summary.tsv                    # Main output table
│   ├── bed/
│   │   ├── srna_0.bed, srna_10.bed, ...
│   │   └── tss_0.bed, tss_10.bed, ...
│   ├── intersects/
│   │   └── 0_0.bed, 0_10.bed, ...     # sRNA ∩ TSS per threshold
│   ├── fasta/
│   │   └── 0_0.fa, 0_10.fa, ...
│   ├── blast/
│   │   └── 0_0.blast, 0_10.blast, ...
│   └── rnie/
│       ├── bed/
│       │   └── 0_0_merged.bed
│       └── fasta/
│           └── 0_0_merged.fa
├── organism2/
│   └── (same structure)
└── ...
```

## Output Interpretation

### Main Output: `summary.tsv`

Tab-separated table with columns:

| Column | Description |
|--------|-------------|
| `score_srna` | Score threshold applied to sRNA predictions |
| `nbr_srna` | Number of sRNA candidates passing threshold |
| `score_tss` | Score threshold applied to TSS predictions |
| `nbr_tss` | Number of TSS passing threshold |
| `intersection_count` | sRNA ∩ TSS (key metric) |
| `blast_hits` | sRNA candidates with BLAST homology |
| `rit_count` | Candidates with predicted RIT |
| `blast_rit` | Highest confidence (sRNA ∩ TSS ∩ RIT + BLAST) |



### System
- Linux/macOS (tested on Ubuntu 20.04, CentOS 7)
- Bash 4.0+
- 32 GB RAM (minimum; 64 GB recommended for multi-organism runs)

### Software
- bedtools ≥ 2.30.0
- BLAST+ ≥ 2.14.1
- SAMtools ≥ 1.15
- Perl ≥ 5.30
- RNIE (with erpin-rho.cm model)
- awk, sed, sort (POSIX standard tools)

## Validation & Testing

The pipeline has been validated on:

**Well-characterized bacteria**:
- *E. coli* K-12 (Enterobacteriales)
- *Salmonella enterica* (Enterobacteriales)
- *Helicobacter pylori* (Epsilonproteobacteria)
- *Staphylococcus aureus* (Bacillales)

**Moderately-characterized bacteria**:
- *Bacillus subtilis* (Bacillales)
- *Mycobacterium tuberculosis* (Actinobacteria)

**Poorly-characterized bacteria**:
- *Methylorubrum extorquens* DM4 (Alphaproteobacteria)
- *Campylobacter jejuni* (Epsilonproteobacteria)
- *Clostridioides difficile* (Firmicutes)


