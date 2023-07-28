
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

## Notes:
There are two accompanying txt files `info/data_db.txt` and `info/human_factor_full_QC.txt` that are important for this pipeline to find the TF and tissue and their bed files.

## To-do:
[ ] - Remove the need for the  `info/data_db.txt` and `info/human_factor_full_QC.txt` file. The user should supply a csv file of the TF, context (tissue) and location of bed files. 