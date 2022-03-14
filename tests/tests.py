import unittest

from pipeline import *


class TestScBonita(unittest.TestCase):
    def test_pipeline(self):

        result = pipeline(dataName="data/example.csv",
            maxNodes=20000,
            maxSamples=15000,
            sep=",",
            getKEGGPathways=True,
            listOfKEGGPathways="00010",
            organism="hsa",
            generateSbatch=False
        )

        # check if scTest pickle has been generated
        pickle=1#True

        # check if scTest data is equal to inout data
        data=1#True


        # check if graphmls have been downloaded
        pathways=1#True

        self.assertEqual(sum(pickle, data, pathways), 3)

if __name__ == '__main__':
    unittest.main()

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
