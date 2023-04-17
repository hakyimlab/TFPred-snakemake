#!/bin/sh
snakemake -s prepare/snakefile.smk --profile profiles/simple/ &&
snakemake -s predict/snakefile.smk -c 1 &&
snakemake -s process/snakefile.smk --profile profiles/cpu/