#!/bin/sh
#SBATCH --partition=debug
#SBATCH -J runSubpops
#SBATCH -o runSubpops.log
#SBATCH -t 1:00:00
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=50G


for d in */ ; do
	echo "$d"
	cd $d
	rm *.graphml
	rm *.log
	rm *.so 
	cd ..
	cp *.sh $d
	cp *.py $d
	cp makefile $d
	cp *.graphml $d
	cp *c $d
	cd $d
	sbatch scBonita_setup.sh
	cd ..
done
