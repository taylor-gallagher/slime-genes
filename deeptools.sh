#!/bin/bash
#SBATCH --job-name=deeptools
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=7-00:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

DEEPTOOLS="/weka/health_sciences/bms/biochemistry/dearden_lab/galta815/containers/deeptools_3.5.3.sif"

for bam in *_sorted.bam; do
    sample_name="${bam%_sorted.bam}"
    
    echo "Processing ${sample_name}..."

    # --binSize 10: High resolution (10bp windows)
    # --normalizeUsing RPKM: Makes different samples comparable
    singularity exec -B /weka "$DEEPTOOLS" bamCoverage \
        -b "$bam" \
        -o "${sample_name}.bw" \
        --binSize 10 \
        --normalizeUsing RPKM \
        --numberOfProcessors 8

done

echo "All BigWig files generated successfully."
