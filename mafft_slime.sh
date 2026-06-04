#!/bin/bash
#SBATCH --job-name=mafft
#SBATCH --cpus-per-task=6
#SBATCH --mem=10G
#SBATCH --time=10:00
#SBATCH --output=mafft.out
#SBATCH --error=mafft.err
#SBATCH --partition=aoraki

singularity exec -B /weka /weka/health_sciences/bms/biochemistry/dearden_lab/galta815/msa/mafft/mafft_7.526.sif \
	mafft --localpair --maxiterate 1000 /weka/health_sciences/bms/biochemistry/dearden_lab/galta815/augustus-ppx/slime_sequences_for_alignment.fasta > slime_sequences_for_augustus.aln
