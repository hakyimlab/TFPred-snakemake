
# TFPred pipeline
## Description: 
This pipeline uses ENFORMER and elastic net to train weights for transcription factor-tissue pairs binding, given the reference genome (fasta file) + intervals to predict on.
## Author: 
Temi
## Date: 
Mon Apr 24 2023
## Usage:
1. Edit the pipeline.yaml file
2. Create your dataframe of TF-tissue pairs --> metadata.txt
3. Activate conda environment
4. Run `bash snakerun.sh` (best to run this in a screen or tmux session)
5. Relax and sip coffee