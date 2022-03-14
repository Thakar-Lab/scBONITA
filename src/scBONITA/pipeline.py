from ruleMaker import *
from keggParser import *
from singleCell import *

# from pathwayAnalysis import *
import glob
import pickle
import argparse
import subprocess


def runAllNetworks(dataFile, maxNodes=20000, maxSamples=50000, partitition='standard', memory = '10G', module = 'anaconda3/2020/07', condaEnv = 'scBonita', pythonVersion = "python3.6", generateSbatch=True):
    for net in glob.glob("*_processed.graphml"):
        if generateSbatch:
            name = net  # [:-8]
            shellHandle = open(name + "_scoreNodes.sh", "w+")
            slurmCommands = str(
                "#!/bin/sh\n#SBATCH --partition=standard\n#SBATCH -J "
                + name
                + "\n#SBATCH -o "
                + name
                + ".log\n#SBATCH -t 24:00:00\n#SBATCH -n 1\n#SBATCH -c 1\n#SBATCH --mem=10G\nmodule load " + str(module) + '\nsource activate ' + str(condaEnv) + '\nmake\n' + str(pythonVersion) + " pipeline.py "
                + "--fullPipeline 0 --dataFile "
                + str(dataFile)
                + " --network "
                + str(net)
                + " --maxNodes "
                + str(maxNodes)
                + " --maxSamples "
                + str(maxSamples)
            )
            shellHandle.write(slurmCommands)
            shellHandle.close()
            shellCommand = name + "_scoreNodes.sh"
            print("Executed ", [shellCommand])
            p = subprocess.Popen(["sbatch", shellCommand])
        else:
            print("For " + net + " run the command: " + str(pythonVerstion) +  " pipeline.py "
                + "--fullPipeline 0 --dataFile "
                + str(dataFile)
                + " --network "
                + str(net)
                + " --maxNodes "
                + str(maxNodes)
                + " --maxSamples "
                + str(maxSamples)
            )

def pipeline(
    dataName="",
    sep=",",
    maxNodes=20000,
    maxSamples=50000,
    getKEGGPathways=True,
    listOfKEGGPathways=[],
    pathwayList=[],
    writeGraphml=True,
    organism="hsa",
    cvThreshold=None,
    binarizeThreshold=0.001,
    generateSbatch=True,
    partitition='standard', 
    memory = '10G',
    module = 'anaconda3/2020/07',
    condaEnv = 'scBonita',
    pythonVersion = 'python3.6'
):
    scTest = singleCell(
        dataName=dataName, sep=sep, maxNodes=maxNodes, maxSamples=maxSamples
    )
    scTest._singleCell__filterData(threshold=cvThreshold)
    if getKEGGPathways:
        pathwayGraphs = find_pathways_kegg(
            geneList=scTest.geneList,
            preDefList=listOfKEGGPathways,
            writeGraphml=writeGraphml,
            organism=organism,
        )
        scTest._singleCell__add_pathways(pathwayGraphs, minOverlap=1)
        print(scTest.pathwayGraphs)
    else:
        scTest._singleCell__add_pathways(pathwayList, minOverlap=1)
    # cut down data size
    nodeIndices = []
    for graph in scTest.pathwayGraphs.keys():
        graph = scTest.pathwayGraphs[graph]
        # collect list of genes from pathways found above
        # find indices of these genes in scObject.geneList
        for node in graph.nodes():
            nodeIndices.append(scTest.geneList.index(node))

    nodeIndices = set(nodeIndices)
    nodeIndices = list(nodeIndices)
    # retain only rows of those genes from binMat
    scTest.expMat = scTest.expMat[nodeIndices, :]
    scTest.binMat = preprocessing.binarize(
        scTest.expMat, threshold=binarizeThreshold, copy=False
    )  # if data >=0: binMat = 1
    scTest.maxNodes = maxNodes
    scTest.maxSamples = maxSamples
    scTest.binMat.resize((scTest.maxNodes, scTest.maxSamples))
    scTest.geneList = [scTest.geneList[node] for node in nodeIndices]
    scTest.nodeList = scTest.geneList
    scTest.nodePositions = [scTest.geneList.index(node) for node in scTest.nodeList]
    pickle.dump(scTest, open(dataName + "scTest.pickle", "wb"))
    runAllNetworks(dataFile=dataName, maxNodes=maxNodes, maxSamples=maxSamples, partition=partition, memory=memory, condaEnv=condaEnv, module=module, pythonVersion=pythonVersion, generateSbatch = generateSbatch)


if __name__ == "__main__":
    # read in arguments from shell scripts
    k = KEGG()
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dataFile",
        help="Specify the name of the file containing processed scRNA-seq data",
        default="",
        type=str,
    )
    parser.add_argument(
        "--fullPipeline",
        type=int,
        choices=[0, 1],
        help="Should scBonita set up the entire pipeline, starting with generation of network topologies?Accepts values 0 or 1.",
    )
    parser.add_argument(
        "--network", help="File name of the network for which rules are to be inferred"
    )
    parser.add_argument(
        "--maxNodes", help="Number of genes in the dataset", default=20000, type=int
    )
    parser.add_argument(
        "--maxSamples", help="Number of cells in the dataset", default=50000, type=int
    )
    parser.add_argument(
        "--separator",
        help="Delimiting character in dataFile. Must be one of , (comma), \s (space) or \t (tab) ",
        default=",",
        type=str,
        choices=[",", "\s", "\t"],
    )
    parser.add_argument(
        "--getKEGGPathways",
        type=int,
        choices=[0, 1],
        help="Should scBonita automatically identify and download KEGG pathways with genes that are in your dataset? You can specify which pathways using the listOfKEGGPathways option, or leave it blank to download all matching KEGG pathways",
        default=0,
        required=False,
    )
    parser.add_argument(
        "--listOfKEGGPathways",
        nargs="+",
        type=str,
        help="Which KEGG pathways should scBonita download? Specify the five letter pathway IDs.",
        required=False,
    )
    parser.add_argument(
        "--pathwayList",
        nargs="+",
        type=str,
        help="Paths to GRAPHML files that should be used for scBONITA analysis. Usually networks from non-KEGG sources, saved in GRAPHML format",
        required=False,
    )
    parser.add_argument(
        "--organism",
        help="Three-letter organism code. Which organism is the dataset derived from?",
        default="hsa",
        type=str,
        choices=k.organismIds,
    )
    parser.add_argument(
        "--cvThreshold",
        help="Minimum coefficient of variation to retain genes for scBONITA analysis",
        default=None,
        type=float,
    )
    parser.add_argument(
        "--binarizeThreshold",
        help="Threshold for binarization of the training dataset. Values above this are classed as 1 (or 'on') and values below this are classed as 0 (or 'off').",
        default=None,
        type=float,
    )
    parser.add_argument(
        "--partition",
        help="SLURM parameter for generated sbatch scripts, if generateSbatch is True",
        default='standard',
        type='str',
    )
    parser.add_argument(
        "--memory",
        help="SLURM parameter for generated sbatch scripts, if generateSbatch is True",
        default='10G',
        type='str',
    )   
    parser.add_argument(
        "--module",
        help="Python/Anaconda module to be loaded in the generated sbatch scripts, if generateSbatch is True",
        default='anaconda3/2020.07',
        type='str',
    )
    parser.add_argument(
        "--condaEnv",
        help="conda environment to be activated in the generated sbatch scripts, if generateSbatch is True",
        default='scBonita',
        type='str',
    )
    parser.add_argument(
        "--pythonVersion",
        help="Python version to be used in the generated sbatch scripts, if generateSbatch is True",
        default='python3.6',
        type='str',
    )  
    results = parser.parse_args()
    fullPipeline = results.fullPipeline
    net = results.network
    maxNodes = results.maxNodes
    maxSamples = results.maxSamples
    dataFile = results.dataFile
    sep = results.separator
    getKEGGPathways = results.getKEGGPathways
    listOfKEGGPathways = results.listOfKEGGPathways
    pathwayList = results.pathwayList
    organism = results.organism
    cvThreshold = results.cvThreshold
    memory = results.memory
    module = results.module
    condaEnv = results.condaEnv
    pythonVersion = results.pythonVersion
    if fullPipeline == 1:
        if dataFile == "":
            dataFile = glob.glob("*.bin")[0]
        else:
            dataFile = str(dataFile)
        print(dataFile)
        pipeline(
            dataName=dataFile,
            maxNodes=maxNodes,
            maxSamples=maxSamples,
            sep=sep,
            getKEGGPathways=getKEGGPathways,
            listOfKEGGPathways=listOfKEGGPathways,
            pathwayList=pathwayList,
            organism=organism,
            cvThreshold=cvThreshold, partitition=partition, memory = memory, module = module, condaEnv = condaEnv, pythonVersion = pythonVersion
        )
    else:
        if fullPipeline == 0:
            if dataFile == "":
                dataFile = glob.glob("*.bin")[0]
            else:
                dataFile = str(dataFile)
            print(dataFile)
            scTest = pickle.load(open(dataFile + "scTest.pickle", "rb"))
            scTest._singleCell__scoreNodes(graph=str(net))

"""
#TEST
# step 1
python3.6 pipeline.py --dataFile "subpop_14.bin" --fullPipeline 1 --maxNodes 20000 --maxSamples 10000 --separator "," --listOfKEGGPathways "00010" "00020" --getKEGGPathways 1 --organism hsa
# step 2
python3.6 pipeline.py --fullPipeline 0 --dataFile subpop_14.bin --network hsa00010.graphml_processed.graphml --maxNodes 20000 --maxSamples 10000
# step 3
python3.6 pipeline.py --fullPipeline 0 --dataFile subpop_14.bin --network hsa00010.graphml_processed.graphml --maxNodes 20000 --maxSamples 10000
# step 4
python3.6 pathwayAnalysis.py --dataFile "subpop_14.bin" --conditions "conditions.txt" --contrast "contrast.txt" --conditions_separator "\t" --contrast_separator "\t"

from pipeline import *
dataName=glob.glob("*.bin")[0]
sep=","
scTest=singleCell(dataName=dataName, sep=sep)
scTest.filterData(threshold=None)
len(scTest.geneList)
find_pathways_kegg(geneList=scTest.geneList, preDefList = ["hsa04670"], writeGraphml=True,  organism="hsa")
scTest.add_pathways(minOverlap=1) #change to use kegg interface
#cut down data size
nodeIndices=[]
for graph in scTest.pathwayGraphs.keys():
    graph = scTest.pathwayGraphs[graph]
    #collect list of genes from pathways found above
    #find indices of these genes in scObject.geneList
    for node in graph.nodes():
        nodeIndices.append(scTest.geneList.index(node))

nodeIndices=set(nodeIndices)
nodeIndices=list(nodeIndices)
#retain only rows of those genes from binMat
print(scTest.expMat.A.shape)
scTest.expMat = scTest.expMat[nodeIndices, :]
print(scTest.expMat.A.shape)
scTest.maxNodes=15000#15000
scTest.maxSamples=10000#60000
scTest.binMat.resize((scTest.maxNodes, scTest.maxSamples)) #restrained to 15000 genes and 10000 samples for testing simulation code, so that I can have a 100 step simulation
scTest.geneList = [scTest.geneList[node] for node in nodeIndices]
scTest.nodeList = scTest.geneList
scTest.nodePositions = [scTest.geneList.index(node) for node in scTest.nodeList]
pickle.dump(scTest, open(dataName+"scTest.pickle", "wb" ))
scTest=pickle.load(open(dataName+"scTest.pickle", "rb"))
for net in glob.glob("*_processed.graphml"):
    name=net[:-8]
    shellHandle=open(name+"_scoreNodes.sh", "w+")
    slurmCommands=str("#!/bin/sh\n#SBATCH --partition=standard\n#SBATCH -J "+name+"\n#SBATCH -o "+name+".log\n#SBATCH -t 24:00:00\n#SBATCH -n 1\n#SBATCH -c 1\n#SBATCH --mem=50G\nmodule load anaconda3/2020.07\nsource activate scBonita\nmake\npython3.6 pipeline.py "+"False "+str(net))
    shellHandle.write(slurmCommands)
    shellHandle.close()
    shellCommand=name+"_scoreNodes.sh"
    print([shellCommand])

graph = glob.glob("*_processed.graphml")[1]
net = graph
scTest.scoreNodes(graph=str(net))

net=nx.read_graphml(graph)
#netGenes = [scTest.geneList.index(gene) for gene in list(net)]
updateBooler=ctypes.cdll.LoadLibrary('./simulator.so')

#Genetic algorithm
population, logbook=scTest.eaMuPlusLambdaAdaptive_altOpt(scSyncBoolC, graph)
out1, out2, model  = scTest.findPopBest(population)
with open(graph+"_rules_GA.txt", "w") as text_file:
    text_file.write(model.writeModel(out2, model))
pickle.dump(out2, open(graph+"_out2.pickle", "wb"))
"""
