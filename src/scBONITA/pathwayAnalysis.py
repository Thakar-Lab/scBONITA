from singleCell import *
from ruleMaker import *
from keggParser import *
from pipeline import *
from attractorAnalysis import *
#import glob
#import seaborn as sns
#import matplotlib.pyplot as plt
#import matplotlib.cm as cm
#from matplotlib.colors import ListedColormap
#from ast import literal_eval
#from random import seed
import argparse

# Perform pathway analysis
#dataFile=glob.glob("*.bin")[0]
#print(dataFile)
#scTest=pickle.load(open(dataFile+"scTest.pickle", "rb"))
#scTest.run_pathway_analysis_from_outs(contrast="contrasts.txt", conditions="conditions.txt", delimiter="\t")

# Perform attractor analysis

#step 0 (quality check): make a plot of equivalent rule set sizes for nodes with an indegree >= 3 and indegree >=2
#tabulate_LS_results() #generates two figures: indegree_3_equivs_histogram.png, indegree_2_equivs_histogram.png

#step 1 - assign cells to network attractors
#assignAttractors(pathwayFiles=None)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataFile", help = "Specify the name of the file containing processed scRNA-seq data", default="", type=str)
    dataFile=results.dataFile
    parser.add_argument("--conditions", help = "Specify the name of the file containing cell-wise condition labels, ie, metadata for each cell in the training dataset. The columns are condition variables and the rows are cells. The first column must be cell names that correspond to the columns in the training data file. The column names must contain the variables specified in the contrast file (see contrast --help for more information).", type=str, required=True)
    conditions=results.conditions
    parser.add_argument("--contrast", help="A text file where each line contains the two conditions (corresponding to column labels in the conditions file) are to be compared during pathway analysis.", required=True, type=str)
    contrast = results.contrast
    parser.add_argument("--conditions_separator", help = "Separator for the conditions file. Must be one of , (comma), \s (space) or \t (tab).", type = str, choices = [",", "\s", "\t"], required=True)
    conditionsSep = results.conditions_separator
    parser.add_argument("--contrast_separator", help = "Separator for the contrast file. Must be one of , (comma), \s (space) or \t (tab).", type = str, choices = [",", "\s", "\t"], required=True)
    contrastSep = results.contrast_separator
    if dataFile == "":
        dataFile=glob.glob("*.bin")[0]
    else:
        dataFile=str(dataFile)
    print(dataFile)
    scTest=pickle.load(open(dataFile+"scTest.pickle", "rb"))
    scTest.run_pathway_analysis_from_outs(contrast=results.contrast, conditions=results.conditions, conditionsSep = conditionsSep, contrastSep = contrastSep)
