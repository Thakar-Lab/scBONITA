#!/bin/sh
#SBATCH --partition=debug
#SBATCH -J setup
#SBATCH -o setup.log
#SBATCH -t 1:00:00
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=60G

module load anaconda3/2020.07

#compile required C code
make

#create a conda environment with all the required packages
#conda env create --name scBonita --file=environment.yml
source activate scBonita

#If the above does not work, try this:
#conda create -n scBonita python=3.6
#conda install -n scBonita scipy scikit-learn numpy seaborn pandas
#conda install -c conda-forge deap
#conda install requests
#conda install -c bioconda bioservices
#conda list --explicit > scBonita_spec_file.txt
#conda create --name scBonita --file scBonita_spec_file.txt


#module load intelpython3/2019.3
#Optionally install required packages
#python -m pip install --user Tkinter matplotlib networkx bioservices deap

#Initiate run
python3.6 reproducible_example_code.py True "blank"
#python3.6 reproducible_example_restartJobs.py
#python3.6 reproducible_example_PA.py

#Check the progress of jobs using:
#squeue -u yourNETID