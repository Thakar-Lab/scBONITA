#!/bin/sh
#SBATCH --partition=debug
#SBATCH -J setup
#SBATCH -o setup.log
#SBATCH -t 1:00:00
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=60G

module load anaconda3/2020.07 #MODIFY THIS

#compile required C code
make

#create a conda environment with all the required packages

source activate scBonita

#Initiate run
# MODIFY condaEnv
python3 pipeline.py --dataFile "Kazer2020expr.csv" --fullPipeline 1 --maxNodes 20000 --separator "," --getKEGGPathways True --organism "hsa" --generateSbatch True --condaEnv "scBonita2" --pythonVersion "python3" --sampleCells True --parallelSearch True --memory "60G" --partition "standard" --binarizeThreshold 0.1

#Check the progress of jobs using:
#squeue -u yourNETID
