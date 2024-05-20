#!/bin/bash

#SBATCH --job-name=train_enpact
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=36:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --error=/project2/haky/temi/projects/lEnpact/logs/train_enpact.err
#SBATCH --output=/project2/haky/temi/projects/lEnpact/logs/train_enpact.log
#SBATCH --account=pi-haky

date
source ~/.bashrc

conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools
export PATH=$PATH:/project2/haky/temi/software/homer/bin

snakemake -s snakefile.smk --configfile config/pipeline.1.yaml --profile profiles/simple/ & \
snakemake -s snakefile.smk --configfile config/pipeline.2.yaml --profile profiles/simple/ & \
snakemake -s snakefile.smk --configfile config/pipeline.3.yaml --profile profiles/simple/ & \
snakemake -s snakefile.smk --configfile config/pipeline.7.yaml --profile profiles/simple/ & \
snakemake -s snakefile.smk --configfile config/pipeline.14.yaml --profile profiles/simple/

