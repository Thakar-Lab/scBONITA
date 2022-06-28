#!/bin/sh
#SBATCH --partition=debug
#SBATCH -J singlePathwaySetup
#SBATCH -o singlePathwaySetup.log
#SBATCH -t 1:00:00
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=60G

module load anaconda3/2020.07

#compile required C code
make

#create a conda environment with all the required packages

source activate scBonita

#Initiate run

python3 pipeline.py --fullPipeline 1 --maxNodes 1000 --separator "," --getKEGGPathways False --pathwayList "temp_series1_net_RELABELED.graphml" --generateSbatch True --partition preempt --time 48:00:00 --binarizeThreshold 0.1 --memory 100G

#python3 pipeline.py --fullPipeline 1 --maxNodes 50000 --separator "," --getKEGGPathways True --listOfKEGGPathways "05110" --generateSbatch True --partition debug --time 1:00:00 --binarizeThreshold 0.001 --memory 60G 

#Check the progress of jobs using:
#squeue -u yourNETID