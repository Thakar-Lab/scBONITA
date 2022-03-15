#!/bin/sh
#SBATCH --partition=debug
#SBATCH -J setup
#SBATCH -o setup.log
#SBATCH -t 1:00:00
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=10G

module load anaconda3/2020.07

#compile required C code
make

#create a conda environment with all the required packages

source activate scBonita


#Initiate run

python3.6 pipeline.py --dataFile "Mild_Severe_Covid_dataset_finalized.csv" --fullPipeline 1 --maxNodes 20000 --maxSamples 15000 --separator "," --listOfKEGGPathways "04066" --getKEGGPathways 1 --organism "hsa"

#Check the progress of jobs using:
#squeue -u yourNETID