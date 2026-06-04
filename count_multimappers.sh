#!/bin/bash
#SBATCH --job-name=count_multimapping
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=7-00:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

BEDTOOLS_SIF="/weka/health_sciences/bms/biochemistry/dearden_lab/galta815/bedtools/bedtools_2.31.1.sif"
BED="clean_loci.bed"

echo "Sample,Tissue,Replicate,Locus,Unique_Reads,Multi_Reads" > mapping_stats_per_locus.csv

for BAM in *_sorted.bam; do
    
    FILENAME=$(basename "$BAM" .bam)
    TISSUE=${FILENAME:0:1}
    REP=${FILENAME:1:1}

    echo "Processing Sample: $FILENAME..."

    while read -r CHROM START END NAME; do
        
        [[ -z "$CHROM" || "$CHROM" == "#"* ]] && continue

        LOCUS_ID="${NAME:-${CHROM}:${START}-${END}}"

        echo -e "$CHROM\t$START\t$END" > current_locus.bed

        
        TOTAL=$(singularity exec --bind /weka "$BEDTOOLS_SIF" bedtools intersect -abam "$BAM" -b current_locus.bed -u | samtools view -c -F 256 -)

        UNIQUE=$(singularity exec --bind /weka "$BEDTOOLS_SIF" bedtools intersect -abam "$BAM" -b current_locus.bed -u | samtools view -F 256 - | grep -c "NH:i:1")

        MULTI=$((TOTAL - UNIQUE))

        echo "$FILENAME,$TISSUE,$REP,$LOCUS_ID,$UNIQUE,$MULTI" >> mapping_stats_per_locus.csv
        
        echo "  - Locus $LOCUS_ID: $TOTAL total reads ($UNIQUE unique)"

    done < "$BED"
done

rm current_locus.bed
echo "Done! Results saved to mapping_stats_per_locus.csv"
