#!/bin/bash
#SBATCH --job-name=generate_profile
#SBATCH --cpus-per-task=36
#SBATCH --mem=100G
#SBATCH --time=1-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=aoraki

perl msa2prfl.pl --qij=default.qij --blockscorefile=slime_gene.blocks.txt slime_sequences_for_augustus.aln > slime_protein.prfl
