
# Enpact pipeline

### Description: 
This pipeline uses ENFORMER and elastic net to train weights for transcription factor binding in a tissue or cell_type, given the reference genome (fasta file) + intervals to predict on. The weights are the effects of the Enformer epigenetic feature on binding status. This method is called Enpact, named after Enformer + IMPACT.

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
There are two accompanying txt files `info/data_db.txt` and `info/human_factor_full_QC.txt` that are important for this pipeline to find the TF and tissue and their peak locations bed files.

## Notebooks:
These contain analysis codes for the pipeline. They are not part of the pipeline itself. You can use these to make diagnostic plots and summaries of the models.

## To-do and Updates
- [ ] Remove the need for the  `info/data_db.txt` and `info/human_factor_full_QC.txt` file. The user should supply a csv file of the TF, context (tissue) and location of bed files. 
- [X] Extend the pipeline to provide summary information of the models, including diagnostic plots.

Tues Aug 22 2023

- [X] Train both linear and logistic models
- [X] Evaluation should be saved into a text file rather than a `.rds` file
- [X] Pipeline can now send jobs to beagle3 (for GPU runs) or caslake as needed
- [X] Added the option to delete ENFORMER predictions on-the-fly as soon as aggregation is done. This will save plenty of storage space when training many models.

Fri Mar 22 2024

- [X] Changes have been made to the way peaks are selected as bound or unbound. Unbound peaks (where there is no motif in a peak) are now selected randomly. For bound peaks:

    1. count the number of bedfiles (experiments) with these peaks
    2. create a probability distribution of the number of experiments that have a peak i.e. larger counts should have higher probabilities of being sub-sampled
    3. select approprately