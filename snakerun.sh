#!/bin/sh
snakemake -s prepare/snakefile.smk --profile profiles/simple/ && #-np &&
snakemake -s predict/snakefile.smk --jobs 10 #--profile profiles/gpu/ &&
# snakemake -s process/snakefile.smk --profile profiles/cpu/