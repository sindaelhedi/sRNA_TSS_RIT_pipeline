 Prerequisites:                                                  
 
✓ sRNA-Detect (outputs: organism_dir/srna.gtf)
✓ TSSAR (outputs: organism_dir/tss.bed)
✓ RNIE (executables: /opt/RNIE/rnie.pl)
✓ RNIE models (files: /opt/RNIE/models/erpin-rho.cm)
✓ BEDTools 
✓ BLAST+ 
✓ SAMtools 
✓ Perl 


# ============================================================
# CONFIG
# ============================================================
BASE_DIR="${1:-}"
OUTPUT_ROOT="${2:-}"


NCPUS="${SLURM_CPUS_PER_TASK:-6}"

TSS_UPSTREAM=60
TSS_DOWNSTREAM=20


srna_scores=(0 10 100 500 1000 1500 2000 2500 3000 3500 10000 20000 40000)
tss_scores=(0 10 100 500 1000 1500 2000 2500 3000 3500 10000 20000 40000)
RNIE_PATH="/path/to/"
CM_MODEL="path/to/"

# ============================================================
# FUNCTIONS
# ============================================================

gtf_to_bed() {
    awk -v s="$3" '$6>=s{print $1,$4,$5,".",$6,$7}' OFS='\t' "$1" \
    | sort -k1,1 -k2,2n > "$2"
}

count_unique() {
    awk '!seen[$1 FS $2 FS $3]++' "$1" | wc -l
}

extend_downstream() {
    awk 'BEGIN{OFS="\t"}{
        name=$1":"$2":"$3":"$6
        if ($6=="+") $3=$3+150
        else $2=($2-150>0)?$2-150:0
        $4=name
        print
    }' "$1" > "$2"
}

run_blast() {
    local fasta="$1"
    local out="$2"
    local db="$3"

    mkdir -p "$(dirname "$out")"
    [[ ! -s "$fasta" ]] && { > "$out"; return; }

    local tmp="${out}.raw"

    # Step 1: full BLAST output (stable)
    blastn -task blastn-short \
        -query "$fasta" \
        -db "$db" \
        -outfmt 6 \
        -num_threads "$NCPUS" \
        > "$tmp" 2>/dev/null || true

    # Step 2: apply filtering AFTER
    awk '!seen1[$1]++ && !seen2[$2]++' "$tmp" > "$out"

    rm -f "$tmp"
}
run_rnie() {
    local fasta="$1"
    local outdir="$2"

    mkdir -p "$outdir"
    [[ ! -s "$fasta" ]] && return

    pushd "$outdir" >/dev/null
    perl "$RNIE_PATH" -m "$CM_MODEL" -f "$fasta" -g --sensitive \
        > rnie.log 2>&1 || true
    popd >/dev/null
}

process_rnie() {
    local gff="$1"
    local out_bed="$2"

    > "$out_bed"
    [[ ! -s "$gff" ]] && return

    awk '$3=="terminator"' "$gff" \
    | awk 'BEGIN{OFS="\t"}{
        split($1,a,"::")
        split(a[1],b,":")
        chr=b[1]; start=b[2]; end=b[3]; strand=b[4]

        split(a[2],c,":|-|\\(|\\)")
        ext_start=c[2]

        term_start=ext_start + $4 - 1
        term_end=ext_start + $5 - 1

        if (strand=="+") {final_start=start; final_end=term_end}
        else {final_start=term_start; final_end=end}

        print chr,final_start,final_end,".",$6,strand
    }' > "$out_bed"
}

# ============================================================
# FULL ANALYSIS
# ============================================================

full_analysis() {

    local tag="$1"
    local bed="$2"
    local genome="$3"
    local db="$4"
    local results="$5"

    mkdir -p "$results"/{fasta,blast,rnie}

    local ext_bed="$results/fasta/${tag}_ext.bed"
    local fasta="$results/fasta/${tag}.fa"

    local rnie_dir="$results/rnie/${tag}"
    local rnie_bed="$results/rnie/${tag}.bed"
    local sorted_bed="$results/rnie/${tag}_sorted.bed"
    local distinct_bed="$results/rnie/${tag}_distinct.bed"
    local distinct_fasta="$results/rnie/${tag}_distinct.fa"

    local blast_out="$results/blast/${tag}.blast"
    local blast_rnie_out="$results/blast/${tag}_rnie.blast"

    # ------------------------
    # FASTA
    # ------------------------
    extend_downstream "$bed" "$ext_bed"

    bedtools getfasta -fi "$genome" -bed "$ext_bed" -s -name \
        -fo "$fasta" 2>/dev/null || true

    run_blast "$fasta" "$blast_out" "$db"

    # ------------------------
    # RNIE
    # ------------------------
    run_rnie "$fasta" "$rnie_dir"

    gff=$(find "$rnie_dir" -name "*rnie.gff" | head -1 || true)
    process_rnie "$gff" "$rnie_bed"

    # ------------------------
    # MERGE RNIE (CRITICAL FIX)
    # ------------------------
    if [[ -s "$rnie_bed" ]]; then

        bedtools sort -i "$rnie_bed" > "$sorted_bed"

        bedtools merge -s -d 2 -i "$sorted_bed" -c 6 -o distinct \
            > "$distinct_bed"

        bedtools getfasta -fi "$genome" -bed "$distinct_bed" \
            -fo "$distinct_fasta" 2>/dev/null || true

        run_blast "$distinct_fasta" "$blast_rnie_out" "$db"
    fi

    # ------------------------
    # COUNTS
    # ------------------------
    blast_n=0
    rit_count=0
    blast_rit=0

    [[ -s "$blast_out" ]] && blast_n=$(wc -l < "$blast_out")
    [[ -s "$distinct_bed" ]] && rit_count=$(wc -l < "$distinct_bed")
    [[ -s "$blast_rnie_out" ]] && blast_rit=$(wc -l < "$blast_rnie_out")

    echo "$blast_n $rit_count $blast_rit"
}

# ============================================================
# MAIN LOOP
# ============================================================

for ORG_DIR in "$BASE_DIR"/*; do

    [[ ! -d "$ORG_DIR" ]] && continue

    ORG=$(basename "$ORG_DIR")
    echo "========== $ORG =========="

    GENOME="$ORG_DIR/genome.fa"
    TSS="$ORG_DIR/tss.bed"
    SRNA="$ORG_DIR/srna.gtf"
    DB_DIR="$ORG_DIR/blast_db"

    BLAST_DB=$(find "$DB_DIR" -name "*.nin" | head -1 | sed 's/.nin$//' || true)

    if [[ ! -f "$GENOME" || ! -f "$TSS" || ! -f "$SRNA" || -z "$BLAST_DB" ]]; then
        echo "[WARNING] Missing files → skipping $ORG"
        continue
    fi

    samtools faidx "$GENOME" 2>/dev/null || true

    RESULTS="$OUTPUT_ROOT/$ORG"
    mkdir -p "$RESULTS"/{bed,intersects,fasta,blast,rnie}

    summary="$RESULTS/summary.tsv"
    echo -e "score_srna\tnbr_srna\tscore_tss\tnbr_tss\tintersection_count\tblast_hits\trit_count\tblast_rit" > "$summary"

    # ------------------------
    # PREP
    # ------------------------
    for s in "${srna_scores[@]}"; do
        gtf_to_bed "$SRNA" "$RESULTS/bed/srna_$s.bed" "$s"
    done

    for t in "${tss_scores[@]}"; do
        awk -v sc="$t" '$5>=sc{print $1,$2,$3,".",$5,$6}' OFS='\t' "$TSS" \
        | sort -k1,1 -k2,2n > "$RESULTS/bed/tss_$t.bed"
    done

    # ------------------------
    # MAIN LOOP
    # ------------------------
    for s in "${srna_scores[@]}"; do

        srna_bed="$RESULTS/bed/srna_$s.bed"
        srna_n=$(count_unique "$srna_bed")
        [[ "$srna_n" -eq 0 ]] && continue

        for t in "${tss_scores[@]}"; do

            tss_bed="$RESULTS/bed/tss_$t.bed"
            tss_n=$(count_unique "$tss_bed")

            intersect="$RESULTS/intersects/${s}_${t}.bed"

            bedtools closest -a "$srna_bed" -b "$tss_bed" -s -D ref \
            | awk -v up="$TSS_UPSTREAM" -v dn="$TSS_DOWNSTREAM" '{
                if ($7 == ".") next
                d=$13
                if ($6=="+" && d>=-up && d<=dn) print $1,$2,$3,$4,$5,$6
                if ($6=="-" && d>=-dn && d<=up) print $1,$2,$3,$4,$5,$6
            }' > "$intersect"

            inter_n=$(count_unique "$intersect")

            blast_n=0
            rit_count=0
            blast_rit=0

            if [[ "$inter_n" -gt 0 ]]; then
                read blast_n rit_count blast_rit <<< "$(full_analysis \
                    "${s}_${t}" "$intersect" "$GENOME" "$BLAST_DB" "$RESULTS")"
            fi

            echo -e "$s\t$srna_n\t$t\t$tss_n\t$inter_n\t$blast_n\t$rit_count\t$blast_rit" >> "$summary"

        done
    done

    echo "[DONE] $ORG"
done

echo "========== ALL DONE =========="




USAGE: bash script.sh <BASE_DIR> <OUTPUT_ROOT>  
ARGUMENTS:
  BASE_DIR      Directory containing organism subdirectories
  OUTPUT_ROOT   Output results directory
  
STRUCTURE:
  BASE_DIR/
  ├── organism1/
  │   ├── genome.fa
  │   ├── tss.bed
  │   ├── srna.gtf
  │   └── blast_db/ (with .nin, .nhr, .nsq files)
  └── organism2/
      ├── genome.fa
      ├── tss.bed
      ├── srna.gtf
      └── blast_db/


RÉSULTATS
 
 
/path/to/results/
├── E_coli/
│   ├── summary.tsv                    ←  PRINCIPAL RESULT
│   ├── bed/
│   │   ├── srna_0.bed, srna_10.bed, ... (sRNA by score)
│   │   └── tss_0.bed, tss_10.bed, ...   (TSS by score)
│   ├── intersects/
│   │   └── 0_0.bed, 0_10.bed, ...      (Intersections sRNA ∩ TSS)
│   ├── fasta/
│   │   └── 0_0.fa, 0_10.fa, ...          (extended sRNA sequences )
│   ├── blast/
│   │   └── 0_0.blast, 0_10.blast, ...  (matched sRNA to databases)
│   └── rnie/
│       ├── 0_0.bed, 0_10.bed, ...      (sRNA+TSS + terminators)
│       └── 0_0_distinct.bed, ...
│
├── M_extorquens/
│   ├── summary.tsv
│   └── (same structure)
│
└── (other organisms...)
 
═══════════════════════════════════════════════════════════════════
 
 summary.tsv
 
 
Colomns :
  score_srna          = score of sRNA
  nbr_srna            = Number of sRNA candidates (sRNA-Detect)
  score_tss           = score of sRNA
  nbr_tss             = Number of sRNA candidates (TSSAR)
  intersection_count  = sRNA ∩ TSS 
  blast_hits          = Validated sRNA by BLAST 
  rit_count           = sRNA ∩ TSS ∩ RIT
  blast_rit           = Validated final sRNA by BLAST


