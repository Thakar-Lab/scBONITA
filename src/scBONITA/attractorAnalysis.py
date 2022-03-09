# import other pieces of scBonita
from singleCell import *
from ruleMaker import *

"""
# import packages
import pickle, re, glob, os, copy, os.path, umap, random, requests, deap, gc
import pandas as pd, seaborn as sns, matplotlib.pyplot as plt, matplotlib.cm as cm, matplotlib.pyplot as plt, scipy.spatial.distance as ssd, networkx as nx, numpy as np, multiprocessing as mp
from os import path
from ast import literal_eval
from random import seed
from time import sleep
from matplotlib.patches import Patch
from scipy.stats import ttest_ind, chi2_contingency
from statsmodels.stats.multitest import multipletests
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

seed(15)  # set seed for reproducibility

# TEST

# Experiment
# Set up
metaData = pd.read_csv("allMetaData.txt", index_col=0)
completeData = pd.read_csv("cell_embeddings.csv", index_col=0)
conditions = pd.read_csv("conditions.txt", sep="\t", index_col=0)
# Reconstruct the singleCell object
objectFile = glob.glob("*.binscTest.pickle")
scObject = pickle.load(open(objectFile[0], "rb"))

# Get network files
networkList = ["hsa0010.graphml_processed.graphml"]

# assign attractors, generate output file
distanceDF = scObject.assignAttractors(
    pathwayFiles=["hsa00010.graphml_processed.graphml"]
)

# Generate UMAP embedding
reducer = umap.UMAP()
embedding = reducer.fit_transform(completeData.values)
plottingData = pd.DataFrame(
    embedding,
    columns=["UMAP dimension 1", "UMAP dimension 2"],
    index=completeData.index,
)
plottingData["Sample"] = [metaData.loc[temp, "batchid"] for temp in plottingData.index]
plottingData["Condition"] = [
    metaData.loc[temp, "PlaqueScore"] for temp in plottingData.index
]

# Make UMAP and attractor frequency plots
scObject.makeAttractorAnalysisPlots(
    plottingData,
    distanceDF,
    allAttractors=False,
    numberOfAttractorsToShow=18,
    cmap="colorblind",
    makeFrequencyPlots=True,
    freqplotsFile="freqPlots_test.pdf",
    makeUMAP=True,
    umapFile="umapPlots_test.pdf",
)
"""
