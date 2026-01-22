#!/bin/bash

# Input arguments
#   1. TSV   : read-to-gene assignment table
#   2. BAM   : original alignment BAM
#   3. EXONS : exon annotation BED
#   4. OUTDIR: output directory
#   5. GENE  : target gene name

TSV=$1
BAM=$2
EXONS=$3
OUTDIR=$4
GENE=$5

mkdir -p "$OUTDIR"

echo "[INFO] Running analysis for $GENE"
echo "[INFO] Output directory: $OUTDIR"

# Step 1: Extract gene-specific paired-end reads

# Extract read1 alignments
samtools view -N <(
    awk -v g="$GENE" '$2==g && /\/1$/{sub(/\/1$/,"",$1); print $1}' "$TSV" | sort -u
) -f 64 "$BAM" -b > "$OUTDIR/read1.bam"

# Extract read2 alignments
samtools view -N <(
    awk -v g="$GENE" '$2==g && /\/2$/{sub(/\/2$/,"",$1); print $1}' "$TSV" | sort -u
) -f 128 "$BAM" -b > "$OUTDIR/read2.bam"

# Merge paired reads
samtools merge -f "$OUTDIR/all_reads.bam" \
    "$OUTDIR/read1.bam" "$OUTDIR/read2.bam"

samtools sort -n "$OUTDIR/all_reads.bam" \
    -o "$OUTDIR/all_reads.sorted.bam"

# Step 2: Define merged genomic region

bedtools bamtobed -i "$OUTDIR/all_reads.sorted.bam" | \
awk '
BEGIN{min=1e18; max=0}
NR==1{chr=$1}
{
    if($2 < min) min=$2;
    if($3 > max) max=$3;
}
END{print chr, min, max}
' OFS='\t' > "$OUTDIR/merged_region.bed"

# Step 3: Intersect merged region with exons

bedtools intersect \
    -a "$EXONS" \
    -b "$OUTDIR/merged_region.bed" \
    > "$OUTDIR/exons_overlap.bed"

# Step 4: Depth calculation

# Exonic depth from original BAM
samtools depth -b "$OUTDIR/exons_overlap.bed" "$BAM" \
    > "$OUTDIR/exon_depth.txt"

# Readthrough depth from gene-specific reads
samtools depth -b "$OUTDIR/merged_region.bed" \
    "$OUTDIR/all_reads.sorted.bam" \
    > "$OUTDIR/readthrough_depth.txt"

# Step 5: Compute split-read ratio (no temp files)

# Identify maximal readthrough depth
max=$(awk 'BEGIN{m=0}{if($3>m)m=$3}END{print m}' \
      "$OUTDIR/readthrough_depth.txt")

# Merge exon depth at positions with maximal readthrough coverage
awk -v m="$max" '
NR==FNR{
    if($3==m) pos[$1":"$2]=$3
    next
}
{
    key=$1":"$2
    if(key in pos) print $0, pos[key]
}' "$OUTDIR/readthrough_depth.txt" "$OUTDIR/exon_depth.txt" \
> "$OUTDIR/merged_depth.txt"

split_ratio=$(awk '{sum+=$4/$3; n++} END{print sum/n}' \
              "$OUTDIR/merged_depth.txt")

# Step 6: Coverage-based readthrough metrics

# Function to compute average depth
avg() { awk '{s+=$3; n++} END{print s/n}' "$1"; }

for region in npepps_remain tbc1d3_readthrough npepps_readthrough tbc1d3_exon1; do
    samtools depth \
        -b "/home/jyguo/tbc/rna_seq/${region}.bed" \
        "$BAM" \
        > "$OUTDIR/${region}_depth.txt"
done

avg1=$(avg "$OUTDIR/npepps_remain_depth.txt")
avg2=$(avg "$OUTDIR/tbc1d3_readthrough_depth.txt")
avg3=$(avg "$OUTDIR/npepps_readthrough_depth.txt")
avg4=$(avg "$OUTDIR/tbc1d3_exon1_depth.txt")

ratioA=$(awk -v a="$avg2" -v b="$avg4" 'BEGIN{print (a-b)/a}')
ratioB=$(awk -v a="$avg1" -v b="$avg3" 'BEGIN{print 1-a/b}')

# Step 7: Output results

cat << EOF > "$OUTDIR/ratios.txt"
Metric	Value
TBC1D3 readthrough_by split-read	$split_ratio
TBC1D3 readthrough_by coverage	$ratioA
NPEPPS readthrough_by coverage	$ratioB
EOF

echo "[DONE] Results saved to: $OUTDIR"