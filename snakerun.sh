#!/bin/sh
snakemake -s prepare/snakefile.smk --configfile config/pipeline.yaml --profile profiles/simple/ -np #&&
# snakemake -s predict/snakefile.smk --configfile config/pipeline.yaml --jobs 10 --profile profiles/gpu/ -np && #--profile profiles/gpu/ &&
# snakemake -s process/snakefile.smk --profile profiles/cpu/

# conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools

# squeue -u $USER | awk '{print $1}' | xargs -n 1 scancel
# snakemake -s predict/snakefile.smk --configfile config/pipeline.yaml --cores 32 --profile profiles/gpu/

# snakemake -s predict/snakefile.smk --configfile config/pipeline.yaml --profile profiles/gpu/


# squeue --format="%.18i %.9P %.30j %.8u %.8T %.10M %.9l %.6D %R" --me

conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools
snakemake -s snakefile.smk --configfile config/pipeline.yaml --profile profiles/simple/ -np