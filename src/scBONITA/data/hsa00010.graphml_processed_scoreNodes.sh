#!/bin/sh
#SBATCH --partition=standard
#SBATCH -J hsa00010.graphml_processed
#SBATCH -o hsa00010.graphml_processed.log
#SBATCH -t 24:00:00
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=10G
module load anaconda3/2020.07
source activate scBonita
make
python3.6 pipeline.py --fullPipeline 0 --dataFile data/trainingData.csv --network hsa00010.graphml_processed.graphml --maxNodes 20000 --maxSamples 10000