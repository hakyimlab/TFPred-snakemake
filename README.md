
# TFPred pipeline

### Description: 
This pipeline uses ENFORMER and elastic net to train weights for transcription factor binding in a tissue or cell_type, given the reference genome (fasta file) + intervals to predict on.

### Author: 
Temi

### Date created: 
Mon Apr 24 2023

## Usage: 
1. Install the conda environment using the file:
    `conda env create -p <<path to env>> -f software/environment.yaml`
2. Edit the `config/pipeline.yaml` file. Instructions are in here.
3. Edit the `config/enformer_base.yaml` file. Instructions are here.
4. Activate conda environment:
    `conda activate <<path to env>>`
5. Run:
    `snakemake -s snakefile.smk --configfile config/pipeline.yaml --profile profiles/simple/`