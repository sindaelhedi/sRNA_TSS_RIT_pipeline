#  sRNA Discovery Pipeline

A comprehensive, generalizable bioinformatic pipeline for systematic discovery and validation of small regulatory RNAs (sRNAs) in bacterial genomes.

---

##  Overview

This pipeline implements a sequential filtering strategy for bacterial sRNA discovery by integrating three independent biological constraints:

1. **sRNA detection** via sRNA-Detect  
2. **Transcription Start Site (TSS) detection** via dRNA-seq (TSSAR)  
3. **Rho-independent terminator (RIT) prediction** via covariance models (RNIE)  

This multi-layer approach significantly reduces transcriptomic noise while enriching high-confidence sRNA candidates.

---

##  Features

- Multi-step biological filtering (TSS + RIT + expression)
- Adaptable to diverse bacterial species
- Parameter tuning for non-model organisms
- Reduces candidate space from thousands to hundreds
- Supports multi-organism analysis
- Structured outputs for downstream validation

---

##  Input

- **Genome sequence** (FASTA)
- **TSS predictions** (BED, from TSSAR)
- **sRNA candidates** (GTF, from sRNA-Detect)
- **BLAST reference database** (nucleotide)

---

##  Output

- **Summary table (TSV)**: candidate metrics across score combinations  
- **Filtered BED files**: sRNA, TSS, intersections  
- **RNIE predictions (GFF)**: intrinsic terminators  
- **BLAST results**: homology validation  
- **FASTA sequences**: for downstream analyses  

---

##  Installation

### Prerequisites

- bedtools ≥ 2.30  
- BLAST+ ≥ 2.14.1  
- SAMtools ≥ 1.15  
- Perl ≥ 5.30  
- RNIE (with covariance models)  

---

### Pipeline usage

```bash
git clone https://github.com/sindalehedi/sRNA-discovery-pipeline.git

##Install sRNA-Detect
git clone https://github.com/BioinformaticsLabAtMUN/sRNA-Detect.git

##TSSAR
wget http://nibiru.tbi.univie.ac.at/TSSAR/download

##Install RNIE 
git clone https://github.com/ppgardne/RNIE.git

### USAGE
```bash
bash pipeline.sh <BASE_DIR> <OUTPUT_ROOT>
