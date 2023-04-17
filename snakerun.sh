
#!/bin/sh

snakemake -s prepare/snakefile.smk --profile profiles/simple/ &&
snakemake -s predict/snakefile.smk -c 1