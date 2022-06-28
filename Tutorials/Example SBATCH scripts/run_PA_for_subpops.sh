#!/bin/sh
#SBATCH --partition=standard
#SBATCH -J runSubpops
#SBATCH -o runSubpops.log
#SBATCH -t 2:00:00
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=10G


for d in */ ; do
    echo "$d"
    cp *.sh $d
    cp *.py $d
    cp makefile $d
    cp *c $d
    #cp kazer_contrast.txt $d
    #cp kazer_metadata.txt $d
    cd $d
    make
    #sbatch reproducible_example_PA.sh
    #sbatch run_pathway_analysis.sh
    make
    module load anaconda3/2020.07
    source activate scBonita
    python3.6 run_pathway_analysis.py
    cd ..
done
