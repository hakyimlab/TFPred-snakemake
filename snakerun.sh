#!/bin/sh
snakemake -s prepare/snakefile.smk --profile profiles/simple/ &&
snakemake -s predict/snakefile.smk -c 1 && #--profile profiles/gpu/
snakemake -s process/snakefile.smk -np -c 1 #--profile profiles/cpu/