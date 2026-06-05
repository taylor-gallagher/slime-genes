#!/bin/bash
#SBATCH --job-name=augustus-ppx
#SBATCH --cpus-per-task=36
#SBATCH --mem=300G
#SBATCH --time=3-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=aoraki

BRAKER_CONFIG="/weka/users/guhjo98p/velvetworm/source_invest/braker_busco_rerun_2026-05-21_existingbam_v3/braker_clean_etp/augustus_config"
SPECIES_NAME="velvetworm_clean_existingbam_v3_20260521" 
IMAGE="/weka/health_sciences/bms/biochemistry/dearden_lab/galta815/augustus-ppx/augustus_v3.3.2dfsg-2-deb_cv1.sif"

singularity exec -B /weka --env AUGUSTUS_CONFIG_PATH=${BRAKER_CONFIG} ${IMAGE} \
        augustus --species=${SPECIES_NAME} \
        --proteinprofile=slime_protein.prfl \
        --predictionStart=745000 \
        --predictionEnd=760000 \
        /weka/health_sciences/bms/biochemistry/dearden_lab/galta815/samtools/Scaffold_18B_Slime_Loci.fasta > Scaffold_18B_Slime_Locus_2.gff
