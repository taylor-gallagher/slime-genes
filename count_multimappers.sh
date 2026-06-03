#!/bin/bash
#SBATCH --job-name=count_multimapping
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=7-00:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

# Define the BED file
BED="slime_gene_loci.bed"

# Print header for the output CSV
echo "Sample,Tissue,Replicate,Unique_Reads,Multi_Reads" > mapping_stats.csv

# Loop through BAM files
for BAM in *_sorted.bam; do
    # Extract Sample Name, Tissue, and Rep from filename 
    SAMPLE=$(basename "$BAM" .bam)
    TISSUE=$(echo $SAMPLE | cut -d'_' -f1)
    REP=$(echo $SAMPLE | cut -d'_' -f2)

    echo "Processing $SAMPLE..."

    samtools view -b -L "$BED" "$BAM" > temp_subset.bam

    # Count Unique Reads (NH=1 and Primary Alignment)
    UNIQUE=$(samtools view -c -F 256 -e 'tag["NH"] == 1' temp_subset.bam)

    #Count Multi-mapping Reads (NH > 1 and Primary Alignment)
    MULTI=$(samtools view -c -F 256 -e 'tag["NH"] > 1' temp_subset.bam)

    # Append to CSV
    echo "$SAMPLE,$TISSUE,$REP,$UNIQUE,$MULTI" >> mapping_stats.csv

    # Clean up temp file
    rm temp_subset.bam
done

echo "Results saved to mapping_stats.csv"
