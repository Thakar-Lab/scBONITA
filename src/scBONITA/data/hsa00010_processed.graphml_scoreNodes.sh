#!/bin/sh
#SBATCH --partition=standard
#SBATCH -J hsa00010_processed.graphml
#SBATCH -o hsa00010_processed.graphml.log
#SBATCH -t 5:00:00
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=10G
module load anaconda3/2020.07
source activate scBonita
make
python3.6 pipeline.py --fullPipeline 0 --dataFile data/trainingData.csv --network hsa00010_processed.graphml --maxNodes 20000 --maxSamples 10000