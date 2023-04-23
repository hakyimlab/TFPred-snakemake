#!/bin/sh
snakemake -s prepare/snakefile.smk --profile profiles/simple/ &&
snakemake -s predict/snakefile.smk --jobs 3 && #--profile profiles/gpu/ &&
snakemake -s process/snakefile.smk --profile profiles/cpu/