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
    #rm *.graphml
    #rm *.log
    #rm *.so 
    cd ..
    #cp *.sh $d
    #cp *.py $d
    #cp makefile $d
    #cp *.graphml $d
    #cp *c $d
    cd $d
    for item in *_scoreNodes.sh
    do
        base="${item%??????????????}"
        out="$base""_importanceScores.csv"
        input="$base""_scoreNodes.sh"
        localSearch="$base""_localErrors1.pickle"
        globalSearch="$base""_out2.pickle"
        [ -f "$input" -a ! -f "$out" ] && echo "$input is missing $out" >&2 #missing importance score calculations
        #[ -f "$input" -a ! -f "$localSearch" ] && echo "$input is missing $localSearch" >&2 #missing local search output
        #[ -f "$input" -a ! -f "$globalSearch" ] && echo "$input is missing $globalSearch" >&2 #missing global search output
        [ -f "$input" -a ! -f "$globalSearch" ] && sbatch $input #(re)execute jobs where global search has not completed
    done
cd ..
done
