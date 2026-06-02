#!/bin/bash
#SBATCH --job-name=tblastn_slime_proteins_transcriptome
#SBATCH --cpus-per-task=36
#SBATCH --mem=200G
#SBATCH --time=7-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=aoraki

tblastn \
	-query slime_protein_sequences.fasta \
	-db /weka/health_sciences/bms/biochemistry/dearden_lab/galta815/blast/trinity-GG-longest-isoforms \
	-out slime_protein_transcriptome_blast_hits.txt \
	-evalue 1e-5 \
	-num_threads 36 \
	-outfmt "6 qseqid sseqid sstart send evalue qcovs length bitscore"
