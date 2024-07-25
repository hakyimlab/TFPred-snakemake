#!/bin/bash

#SBATCH 
#SBATCH --partition=beagle3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=enformer_predict
#SBATCH --account=pi-haky
#SBATCH --time=06:00:00
#SBATCH --gres=gpu:2

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools

nvidia-smi

# make sure that in your config, provider is set to "local" before running this script
python ${1} --parameters ${2}