#!/bin/bash
#SBATCH --job-name=fastBlockSearch
#SBATCH --cpus-per-task=36
#SBATCH --mem=300G
#SBATCH --time=1-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=aoraki

singularity exec -B /weka /weka/health_sciences/bms/biochemistry/dearden_lab/galta815/augustus-ppx/augustus_v3.3.2dfsg-2-deb_cv1.sif \
	fastBlockSearch --cutoff=1.1 /weka/health_sciences/bms/biochemistry/dearden_lab/galta815/augustus-ppx/Slime_Loci.fasta slime_protein.prfl \
	2> fastblock.out

