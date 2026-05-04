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
CM_MODEL="/path/to/"

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

# ============================================================
# CONVERT SPACES TO TABS IN BED FILE
# ============================================================
# bedtools requires TAB-separated files, not spaces

fix_bed_format() {
    local bed="$1"
    local output="$2"
    
    # Convert multiple spaces/whitespace to tabs
    sed 's/[[:space:]]\+/\t/g' "$bed" > "$output"
}

run_blast() {
    local fasta="$1"
    local out="$2"
    local db="$3"

    mkdir -p "$(dirname "$out")"
    [[ ! -s "$fasta" ]] && { > "$out"; return; }

    local tmp="${out}.raw"
    local seq_count=$(grep -c '^>' "$fasta" 2>/dev/null || echo 0)
    echo "[BLAST] Input: $seq_count sequences from $(basename "$fasta")" >&2

    blastn -task blastn-short \
        -query "$fasta" \
        -db "$db" \
        -outfmt 6 \
        -num_threads "$NCPUS" \
        -evalue 1e-1 \
        -perc_identity 95 \
        > "$tmp" 2>/dev/null || true

    local raw_count=$(wc -l < "$tmp" 2>/dev/null || echo 0)
    echo "[BLAST] Raw hits: $raw_count" >&2

    # CORRECT FILTERING: Keep only first hit per query sequence
    awk '!seen[$1]++' "$tmp" > "$out"
    
    local filtered_count=$(wc -l < "$out" 2>/dev/null || echo 0)
    echo "[BLAST] After filtering: $filtered_count" >&2

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
        chr=b[1]
        start=b[2]
        end=b[3]
        strand=b[4]

        split(a[2],c,":|-|\\(|\\)")
        ext_start=c[2]

        term_start=ext_start + $4 - 1
        term_end=ext_start + $5 - 1

        if (strand=="+") {
            final_start=start
            final_end=term_end
        } else {
            final_start=term_start
            final_end=end
        }

        key = chr ":" start ":" end ":" strand

        if (!(key in best_score) || $6 > best_score[key]) {
            best_score[key] = $6
            best_data[key] = chr "\t" final_start "\t" final_end "\t.\t" $6 "\t" $7
        }
    }
    END {
        for (k in best_data) {
            print best_data[k]
        }
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

    local bed_count=$(wc -l < "$bed")
    echo "[ANALYSIS] Processing ${tag}: ${bed_count} intersected sRNAs" >&2

    # ========================
    # FIX BED FORMAT (spaces → tabs)
    # ========================
    
    local bed_fixed="$results/fasta/${tag}_fixed.bed"
    fix_bed_format "$bed" "$bed_fixed"
    
    local fixed_count=$(wc -l < "$bed_fixed")
    echo "[ANALYSIS] Fixed BED format: ${fixed_count} lines" >&2

    # ========================
    # PART 1: BLAST ANALYSIS
    # ========================
    # BLAST should use ORIGINAL sRNA sequences (not extended)
    
    local fasta="$results/fasta/${tag}.fa"
    
    # Extract FASTA directly from FIXED intersected BED
    bedtools getfasta -fi "$genome" -bed "$bed_fixed" -s -name -fo "$fasta" 2>/dev/null || true

    local extracted_count=$(grep -c '^>' "$fasta" 2>/dev/null || echo 0)
    echo "[ANALYSIS] Extracted original FASTA: ${extracted_count} sequences" >&2

    # BLAST the original sRNA sequences
    local blast_out="$results/blast/${tag}.blast"
    run_blast "$fasta" "$blast_out" "$db"

    # ========================
    # PART 2: RNIE ANALYSIS
    # ========================
    # RNIE needs to find TERMINATORS downstream
    
    local ext_bed="$results/fasta/${tag}_ext.bed"
    local ext_fasta="$results/fasta/${tag}_ext.fa"
    
    # Extend coordinates for terminator search
    extend_downstream "$bed_fixed" "$ext_bed"
    
    # Extract extended FASTA (includes downstream region)
    bedtools getfasta -fi "$genome" -bed "$ext_bed" -s -name -fo "$ext_fasta" 2>/dev/null || true

    local ext_extracted=$(grep -c '^>' "$ext_fasta" 2>/dev/null || echo 0)
    echo "[ANALYSIS] Extracted extended FASTA: ${ext_extracted} sequences" >&2

    # Run RNIE on extended sequences to find terminators
    local rnie_dir="$results/rnie/${tag}"
    run_rnie "$ext_fasta" "$rnie_dir"

    local rnie_bed="$results/rnie/${tag}.bed"
    local gff=$(find "$rnie_dir" -name "*rnie.gff" | head -1 || true)
    process_rnie "$gff" "$rnie_bed"

    # ========================
    # PART 3: RNIE VALIDATION
    # ========================
    
    local sorted_bed="$results/rnie/${tag}_sorted.bed"
    local distinct_bed="$results/rnie/${tag}_distinct.bed"
    local distinct_fasta="$results/rnie/${tag}_distinct.fa"
    local blast_rnie_out="$results/blast/${tag}_rnie.blast"

    if [[ -s "$rnie_bed" ]]; then
        bedtools sort -i "$rnie_bed" > "$sorted_bed"
        bedtools merge -s -d 2 -i "$sorted_bed" -c 6 -o distinct > "$distinct_bed"
        bedtools getfasta -fi "$genome" -bed "$distinct_bed" -fo "$distinct_fasta" 2>/dev/null || true
        
        # BLAST the RNIE-derived sequences
        run_blast "$distinct_fasta" "$blast_rnie_out" "$db"
    fi

    # ========================
    # RETURN COUNTS
    # ========================
    local blast_n=0
    local rit_count=0
    local blast_rit=0

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

    # ========================
    # PREP: Create BED files
    # ========================
    for s in "${srna_scores[@]}"; do
        gtf_to_bed "$SRNA" "$RESULTS/bed/srna_$s.bed" "$s"
    done

    for t in "${tss_scores[@]}"; do
        awk -v sc="$t" '$5>=sc{print $1,$2,$3,".",$5,$6}' OFS='\t' "$TSS" \
        | sort -k1,1 -k2,2n > "$RESULTS/bed/tss_$t.bed"
    done

    # ========================
    # MAIN SCORING LOOP
    # ========================
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





