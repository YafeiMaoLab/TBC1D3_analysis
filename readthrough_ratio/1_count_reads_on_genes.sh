#!/bin/bash

# This script is used to count the number of reads on each gene in a BAM file.

# input:
# 1. GFF file
# 2. BAM file
# 3. output prefix

# output:
# 1. TSV file with the number of reads on each gene

# Description:
# - unmapped: unmapped reads
# - unspcified: aligned to non-gene regions (introns, intergenic, etc.)
# - [GeneName-...-GeneName]: assigned to multiple genes:
# - [GeneName]: assigned to unique genes

gff=$1
bam=$2
outprefix=$3

module load samtools
module load gffread

# Exit immediately on error
set -euo pipefail

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <annotation.gff> <input.bam> <output_prefix>"
    exit 1
fi

# index BAM
if [ ! -f "${bam}.bai" ]; then
    samtools index "$bam"
fi

# make directory for output if not exists
outdir=$(dirname "$outprefix")  # get the parent directory
rm -rf "$outdir"
mkdir -p "$outdir" 

# sort bam
samtools sort -o "${outprefix}.sorted.bam" "$bam"
samtools index "${outprefix}.sorted.bam"

# extract mapped reads
samtools view -F 260 -b "${outprefix}.sorted.bam" > "${outprefix}.mapped.bam" # -F 260
samtools view -F 4 -b "${outprefix}.sorted.bam" > "${outprefix}.unmapped.bam"
samtools fastq "${outprefix}.unmapped.bam" | gzip -c > "${outprefix}.unmapped.fastq.gz"
# -F 4: exclude unmapped reads
# -F 260: exclude secondary alignments and unmapped reads
# -F 308: primary only, no supplementary

# count total reads
echo "[INFO] Counting total reads..."
total=$(samtools view -c -F 2304 "${outprefix}.sorted.bam") # only include primary alignments and unaligned reads
mapped=$(samtools view -c -F 2308 "${outprefix}.sorted.bam")
unmapped=$((total - mapped))

# get gene exon regions from GFF into bed format
echo "[INFO] Extracting gene exon regions from GFF..."
gffread --bed "$gff" > "${outprefix}.gene.bed"
# for each gene, get the exon regions
awk -F'\t' -v OFS='\t' -v outdir="$outdir" '
{
    chr=$1
    strand=$6
    exonCount=$10
    exonLengths=$11
    exonStarts=$12
    geneID="unknown"
    # parse geneID from last column
    if($13 ~ /geneID=/){
        match($13,/geneID=([^;]+)/,a)
        if(a[1]!="") geneID=a[1]
    }

    split(exonLengths, len_arr, ",")
    split(exonStarts, start_arr, ",")

    for(i=1;i<=exonCount;i++){
        if(len_arr[i]!=""){
            start = $7 + start_arr[i]       # chrom start = chromStart + exonStart
            end   = start + len_arr[i]      # chrom end = chrom start + exonLength
            exon_name = geneID "-" $4 "-exon" i
            print chr, start, end, exon_name >> outdir"/"geneID".exons.bed"
        }
    }
}' "${outprefix}.gene.bed"

# generate merged bed file
if [ -f "${outdir}/merged.exons.bed" ]; then
    rm "${outdir}/merged.exons.bed"
fi
cat "${outdir}"/*.exons.bed > "${outdir}/merged.exons.bed"

# assign reads to genes
echo "[INFO] Assigning reads to genes..."
bedtools intersect -a ${outprefix}.mapped.bam -b ${outdir}/merged.exons.bed -v > ${outprefix}.unspecified.bam # filter unspecified read
###bedtools intersect -a ${outprefix}.mapped.bam -b ${outdir}/merged.exons.bed -wa | \
###samtools view -h - | awk '$1 ~ /^@/ || !seen[$1]++' | samtools view -b -o ${outprefix}.assigned.bam # count assigned reads (each read appears only once)
bedtools intersect -a ${outprefix}.mapped.bam -b ${outdir}/merged.exons.bed -wa -u > ${outprefix}.assigned.bam
unspecified=$(samtools view -c -F 2308 ${outprefix}.unspecified.bam)
assigned=$(samtools view -c -F 2308 ${outprefix}.assigned.bam)

# Get overlaps of reads with exons
bedtools intersect -a "${outprefix}.assigned.bam" -b "${outdir}/merged.exons.bed" -wa -wb -bed > "${outprefix}.read_overlaps.tsv"

# Assign reads to unique/multiple genes
awk -F'\t' '
{
    read=$4
    split($16, arr, "-")          # exon_name -> GeneName-GeneID-exonN
    gene=arr[1]
    reads[read][gene]=1
}
END{
    for(r in reads){
        n = 0
        genes = ""
        for(g in reads[r]){
            n++
            genes=(genes=="" ? g : genes"-"g)
        }
        # MODIFIED: each read assigned to exactly one label (single-gene or multi-gene)
        print r "\t" genes
    }
}' "${outprefix}.read_overlaps.tsv" > "${outprefix}.read_assignment.tsv"

# Assigned
mkdir -p "${outdir}/assigned"
###while read -r read name; do
###    samtools view -b -h "${outprefix}.mapped.bam" | \
###        awk -v r="$read" 'BEGIN{OFS="\t"} $1==r || /^@/' >> "${outdir}/assigned/${name}.bam"
###done < "${outprefix}.read_assignment.tsv"

awk -F'\t' -v OFS='\t' '{print $2, $1}' "${outprefix}.read_assignment.tsv" | sort -k1,1 > "${outdir}/gene_read_pairs.tsv"
samtools view -F 2308 ${outprefix}.assigned.bam | awk '{print $1}' | sort -u > ${outdir}/valid_read_names.txt
awk 'NR==FNR {keep[$1]; next} ($2 in keep)' ${outdir}/valid_read_names.txt "${outdir}/gene_read_pairs.tsv" > "${outdir}/gene_read_pairs.filtered.tsv"

cut -f1 "${outdir}/gene_read_pairs.filtered.tsv" | sort | uniq | while read gene; do
    safe_name=$(echo "$gene" | tr '/<>:"|?* ' '_')   # 安全文件名
    grep -P "^$gene\t" "${outdir}/gene_read_pairs.filtered.tsv" | cut -f2 > "${outdir}/assigned/${safe_name}.readlist.txt"
    samtools view -b -h -N "${outdir}/assigned/${safe_name}.readlist.txt" "${outprefix}.mapped.bam" -o "${outdir}/assigned/${safe_name}.bam"
    rm "${outdir}/assigned/${safe_name}.readlist.txt"
done

# Output to summary file
touch "${outprefix}.summary.txt"
echo -e "--------------------------------" >> "${outprefix}.summary.txt"
echo -e "# Total = Mapped + Unmapped" >> "${outprefix}.summary.txt"
echo -e "Total\t${total}" >> "${outprefix}.summary.txt"
echo -e "Unmapped\t${unmapped}" >> "${outprefix}.summary.txt"
echo -e "Mapped\t${mapped}" >> "${outprefix}.summary.txt"
echo -e "--------------------------------" >> "${outprefix}.summary.txt"
echo -e "# Mapped = Unspecified + Assigned" >> "${outprefix}.summary.txt"
echo -e "Unspecified\t${unspecified}" >> "${outprefix}.summary.txt"
echo -e "Assigned\t${assigned}" >> "${outprefix}.summary.txt"
echo -e "--------------------------------" >> "${outprefix}.summary.txt"
echo -e "# Assigned = Sum of assigned reads on each label" >> "${outprefix}.summary.txt"
for bam in `ls ${outdir}/assigned/*.bam`; do
    name=$(basename $bam .bam)
    num=$(samtools view -c -F 2308 $bam)
    echo -e "${name}\t${num}" >> "${outprefix}.summary.txt"
done

# collect bam
mkdir -p "${outdir}/bam"
mv ${outdir}/*.bam "${outdir}/bam"
mv ${outdir}/assigned/*.bam "${outdir}/bam"
rm -r ${outdir}/assigned
