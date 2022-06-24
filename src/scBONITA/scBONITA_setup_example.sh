#!/bin/sh
#SBATCH --partition=debug
#SBATCH -J setup
#SBATCH -o setup.log
#SBATCH -t 0:30:00
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=10G

module load anaconda3/2020.07

#compile required C code
make

#create a conda environment with all the required packages

source activate scBonita


#Initiate run

python3.6 pipeline.py --dataFile "mildsevereCOVID__normalized_integrated_counts.csv" --fullPipeline 1 --maxNodes 20000 --maxSamples 100 --separator "," --listOfKEGGPathways "04066" "04010" "04150" "04020" "00020" "04370" "04630" "04668" "04060" "04514" "04670" "04512" "00010" "04625" "04662" "04120" "04062" --getKEGGPathways True --organism "hsa" --generateSbatch True

#Check the progress of jobs using:
#squeue -u yourNETID