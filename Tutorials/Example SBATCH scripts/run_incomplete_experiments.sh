#!/bin/sh
#SBATCH --partition=debug
#SBATCH -J run_incompletes
#SBATCH -o run_incompletes.log
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
    for item in *_scoreNodes.sh
    do
        base="${item%??????????????}"
        out="$base""_importanceScores.csv"
        input="$base""_scoreNodes.sh"
        [ -f "$input" -a ! -f "$out" ] && echo "$input is missing $out" >&2
        sbatch $input
    done
cd ..
done
