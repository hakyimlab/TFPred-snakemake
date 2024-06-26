
# Enpact pipeline

### Description: 
This pipeline uses ENFORMER and elastic net to train weights for transcription factor binding in a tissue or cell_type, given the reference genome (fasta file) + intervals to predict on. The weights are the effects of the Enformer epigenetic feature on binding status. This method is called Enpact, named after Enformer + IMPACT.

### Author: 
Temi

### Date created: 
Mon Apr 24 2023

## Usage: 

### Notes:

There are 2 ways to use this pipeline. 

1. Train Enpact models using a reference epigenome (from Enformer): Here, the pipeline will not run Enformer and will instead make use of the reference epigenome provided. This is faster than the second approach below. To use this approach, you need to provide the reference epigenome in the [config/pipeline.yaml](./minimal/pipeline.minimal.yaml) file, and set the parameter `run_enformer` to `false`.

2. Train Enpact models using Enformer: Here, the pipeline will run Enformer to generate the reference epigenome. This is slower than the first approach above. To use this approach, you need to set the parameter `run_enformer` to `true` in the [config/pipeline.yaml](./minimal/pipeline.minimal.yaml) file.


### Helpers:

There is a notebook [here](./notebooks/prepare_samples.qmd) to help generate the inputs. 

Otherwise, you can look at minimal examples of the [models.data.yaml](./minimal/models.data.yaml) and [models.run.tsv](./minimal/models.run.tsv) files to recreate yours.

### Input:

1. Install the conda environment using the file:
    `conda env create -p <<path to env>> -f software/environment.yaml`
2. Edit the `config/pipeline.yaml` file. Instructions are in here.
3. Edit the `config/enformer_base.yaml` file. Instructions are here.
4. Activate conda environment:
    `conda activate <<path to env>>`
5. Run:
    `snakemake -s snakefile.smk --configfile config/pipeline.yaml --profile profiles/simple/`

## Notebooks:
These contain analysis codes for the pipeline. They are not part of the pipeline itself. You can use these to make diagnostic plots and summaries of the models.

## To-do and Updates
- [X] Remove the need for the  `info/data_db.txt` and `info/human_factor_full_QC.txt` file. The user should supply a csv file of the TF, context (tissue) and location of bed files. 
- [X] Extend the pipeline to provide summary information of the models, including diagnostic plots.

Tues Mar 27 2024

- [X] Added the option to test models on an held-out chromosome. This is the default option. Set to false if you want to use random motifs across the genome.

Tues Mar 26 2024

- [X] Intersection of peaks is now calculated using `bedtools intersect`. This is faster than using code in R that calculates the intersection of peaks.

Fri Mar 22 2024

- [X] Changes have been made to the way peaks are selected as bound or unbound. Unbound peaks (where there is no motif in a peak) are now selected randomly. For bound peaks:

    1. count the number of bedfiles (experiments) with these peaks i.e. binding counts
    2. assign probabilities to each unique binding count that have a peak i.e. larger counts should have higher probabilities of being sub-sampled, and lower counts should have lower probabilities
    3. select approprately

Tues Aug 22 2023

- [X] Train both linear and logistic models
- [X] Evaluation should be saved into a text file rather than a `.rds` file
- [X] Pipeline can now send jobs to beagle3 (for GPU runs) or caslake as needed
- [X] Added the option to delete ENFORMER predictions on-the-fly as soon as aggregation is done. This will save plenty of storage space when training many models.

Sun Jun 9 2024

- [X] Extensive modification to how the pipeline should run.
