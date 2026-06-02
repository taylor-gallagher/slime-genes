#!/bin/bash
#SBATCH --job-name=samtools_index
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=7-00:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

for bam in *_sorted.bam; do
    echo "Indexing $bam with CSI..."
    # -c forces CSI format, -@ 4 uses 4 threads
    samtools index -c -@ 4 "$bam"
done
