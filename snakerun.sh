#!/bin/sh
snakemake -s snakefile.smk --configfile config/pipeline.yaml --profile profiles/simple/ -np #&&
# snakemake -s predict/snakefile.smk --configfile config/pipeline.yaml --jobs 10 --profile profiles/gpu/ -np && #--profile profiles/gpu/ &&
# snakemake -s process/snakefile.smk --profile profiles/cpu/

# conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools

# squeue -u $USER | awk '{print $1}' | xargs -n 1 scancel
# snakemake -s predict/snakefile.smk --configfile config/pipeline.yaml --cores 32 --profile profiles/gpu/

# snakemake -s predict/snakefile.smk --configfile config/pipeline.yaml --profile profiles/gpu/


# squeue --format="%.18i %.9P %.40j %.8u %.8T %.10M %.9l %.6D %R" --me
screen
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools
export PATH=$PATH:/beagle3/haky/users/temi/software/homer/bin

# export PATH=$PATH:/beagle3/haky/users/temi/software/homer/bin

snakemake -s snakefile.smk --configfile config/pipeline.yaml --profile profiles/simple/ --resources load=50 -F >> run.txt 2>&1  


snakemake -s snakefile.smk --configfile config/pipeline.yaml --profile profiles/simple/ --resources load=50 --report reports/report.html -F >> run.txt 2>&1  

# Generate report
snakemake -s snakefile.smk --configfile config/pipeline_2023-12-01_enformer_vs_tfpred.yaml --profile profiles/simple/ --resources load=50 -F >> run_2023-12-01_enformer_vs_tfpred.txt 2>&1

# snakemake -s snakefile.smk --configfile config/pipeline.yaml --profile profiles/simple/ -R prepare_training_data

# snakemake -s snakefile.smk --configfile config/pipeline.yaml --profile profiles/simple/ --resources load=50 --report reports/report.html -F --rerun-incomplete >> run.txt 2>&1  
# snakemake -s snakefile.smk --configfile config/pipeline.yaml --profile profiles/simple/ --resources load=50 --report reports/report.html >> run.txt 2>&1  

snakemake -s snakefile.smk --configfile config/pipeline.yaml --profile profiles/simple/ -R create_training_set -np

snakemake -s snakefile.smk --configfile config/pipeline.yaml --profile profiles/simple/ -R create_training_set -np

snakemake -s snakefile.smk --configfile config/pipeline_2023-12-01_enformer_vs_tfpred.yaml --profile profiles/simple/ --resources load=50 >> run_2023-12-01_enformer_vs_tfpred.txt 2>&1

snakemake -s snakefile.smk --configfile config/pipeline.yaml --profile profiles/simple/ -np



snakemake -s snakefile.smk --configfile config/pipeline.enpact.yaml --profile profiles/simple/ -np > dryrun.out


snakemake -s snakefile.smk --configfile experiments/bpnet/pipeline.bpnet.yaml --profile profiles/simple/ -np > dryrun.out


snakemake -s snakefile.smk --configfile config/pipeline.734.yaml --profile profiles/simple/ -np > dryrun.out