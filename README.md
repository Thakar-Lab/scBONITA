![GitHub last commit](https://img.shields.io/github/last-commit/Thakar-Lab/scBONITA?style=for-the-badge)

# scBonita: single-cell Boolean Omics Network Invariant Time Analysis

## Publication: 

Medrxiv link coming soon

## Description

To use single-cellRNA seq data to develop executable models of signaling pathways that drive cellular states, we developed the single-cell Boolean Omics Network Invariant Time Analysis (scBONITA) algorithm. ScBONITA 
* uses scRNA-seq data to infer Boolean regulatory rules for topologically characterized networks 
* prioritizes genes based on their impact on signaling
* performs pathway analysis, and 
* maps sequenced cells to characteristic signaling states of these networks

We show in our manuscript that scBONITA learns dynamic models of signaling pathways to identify pathways significantly dysregulated in HIV-associated atherosclerosis. These pathways indicate that cell migration into the vascular endothelium is significantly affected by dysregulated lipid signaling in AS+ PLWH. Dynamic modeling facilitates pathway-based characterization of cellular states that are not apparent in gene expression analyses.

## Graphical abstract

![alt text](https://github.com/Thakar-Lab/scBONITA/blob/main/graphical_abstract_scBonita.png?raw=true)

## Keywords

single-cell RNA sequencing; Boolean networks; HIV; atherosclerosis; pathway analysis; network modeling

## Contact: 
Juilee_Thakar@URMC.Rochester.edu; Mukta_Palshikar@URMC.Rochester.edu

## Bug Tracker:

https://github.com/Thakar-Lab/scBONITA/issues

## Installation instructions

**scBONITA requires Python 3.6 or newer.**

scBONITA is currently designed to be run on SLURM systems, with minimally complicated transfer to other distributed computing systems. 

We recommend that the scBONITA pipeline is used on a high-performance system with adequate computational resources. We find that the amount of RAM required for large scRNA-seq datasets is in excess of 10G; we hence don't recommend running rule inference and attractor analysis on a desktop computer.


### Use package from GitHub

1. Download a zipped folder from github.com/Thakar-Lab/scBONITA, OR

1. Use git to clone this repository
    `git clone https://github.com/YOUR-USERNAME/YOUR-REPOSITORY`

Modify the C file to reflect the size of the training data set
    - line 7: NODE is the number of nodes in the dataset (or larger)
    - line 8: CELL is the number of cells in the dataset (or larger)

Next, the C code must be compiled using the make file. Navigate to the folder in which the scBONITA package is located and in a Bash terminal, type
    `make`

Install the scBonita conda environment from the provided yml file, and activate the conda environment.

## Step 1: Set up the scBONITA pipeline.

scBONITA needs a training dataset in matrix-style forma; this is usually a tab or comma-delimited file with columns as cells and rows as features. The first column should be feature names and the first row should be cell IDs. The units of the expression data will typically be a variant of log2(TPM +1).


The setup pipeline has the following parameters:
1. dataFile: Specify the name of the file containing processed scRNA-seq data
1. fullPipeline Should scBonita set up the entire pipeline, starting with generation of network topologies? Accepts values 0 or 1
1. network File name of the network for which rules are to be inferred
1. maxNodes Number of genes in the dataset
1. maxSamples Number of cells in the dataset
1. separator Delimiting character in dataFile. Must be one of , (comma), \s (space) or \t (tab)
1. getKEGGPathways Should scBonita automatically identify and download KEGG pathways with genes that are in your dataset? You can specify which pathways using the 
1. listOfKEGGPathways option, or leave it blank to download all matching KEGG pathways
1. listOfKEGGPathways Which KEGG pathways should scBonita download? Specify the five letter pathway IDs.
1. organism Three-letter organism code. Which organism is the dataset derived from?
1. cvThreshold: Minimum coefficient of variation to retain genes for scBONITA analysis
The following line of code runs scBONITA setup for a 20000*10000 comma-separated data set "example.csv". It downloads the KEGG pathways hsa00010 and hsa00020 for rule inference.

`python3.6 pipeline.py --dataFile "example.csv" --fullPipeline 1 --maxNodes 20000 --maxSamples 10000 --separator "," --listOfKEGGPathways "00010" "00020" --getKEGGPathways 1 --organism hsa cvThreshold None`

### Step 2: Rule inference and calculation of node importance score for the networks specified in Step 1 (setup).

Step 1 generates sbatch files that enter the specified slurm queue. In a typical use case, these jobs should execute automatically. We recommend that users periodically check the slurm queue and the log files.

### Step 3: Pathway Analysis

To perform pathway analysis, scBONITA uses the rules generated in Step 2. In addition, scBONITA requires (a) a metadata file specifiying the treatments/experimental variables for each cell and (b) a contrast file specifying the pairs of treatments to be compared.

The pathway analysis script has the following arguments:

1. dataFile Specify the name of the file containing processed scRNA-seq data
1. conditions Specify the name of the file containing cell-wise condition labels, ie, metadata for each cell in the training dataset. The columns are condition variables and the rows are cells. The first column must be cell names that correspond to the columns in the training data file. The column names must contain the variables specified in the contrast file (see contrast --help for more information).
1. contrast A text file where each line contains the two conditions (corresponding to column labels in the conditions file) are to be compared during pathway analysis.
1. conditions_separator Separator for the conditions file. Must be one of , (comma), \s (space) or \t (tab).
1. contrast_separator Separator for the contrast file. Must be one of , (comma), \s (space) or \t (tab).

For example:

`python3.6 pathwayAnalysis.py --dataFile "example.csv" --conditions "conditions.txt" --contrast "contrast.txt" --conditions_separator "\t" --contrast_separator "\t"`

** Please refer to the tutorial "Pathway_Analysis_With_scBONITA.ipynb" for suggestions on analysis of the output of scBONITA PA **

### Step 4: Attractor Analysis

** Please refer to the tutorial "Attractor_Analysis_With_scBONITA.ipynb" for a demonstration of attractor analysis with scBONITA using a test dataset packaged with scBONITA** 

