#!/bin/sh
#SBATCH --partition=preempt
#SBATCH -J PA
#SBATCH -o pathwayAnalysis.log
#SBATCH -t 0:20:00
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=10G
module load anaconda3/2020.07
source activate scBonita
make
python run_pathway_analysis.py