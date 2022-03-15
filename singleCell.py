from ruleMaker import *
from keggParser import *
import pickle
import scipy.sparse as sparse
from scipy.stats.stats import spearmanr
import glob
import numpy as np
from os import path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from ast import literal_eval
import random
import os
import scipy
from statsmodels.stats.multitest import multipletests
from pathlib import Path


class singleCell(ruleMaker):

    """Class for single-cell experiments"""

    def __init__(
        self, dataName, sep, maxNodes=15000, maxSamples=10000, binarizeThreshold=0.001, sampleCells=True
    ):
        """Read in pre-processed data and binarize by threshold"""
        data = np.loadtxt(dataName, delimiter=sep, dtype="str")

        if maxSamples > 15000 or sampleCells:
            sampledCellIndices = self.__sampleCells(data=data, number_cells = max(15000, maxSamples))
            self.geneList = data[1:, 0]
            self.sampleList = data[0, sampledCellIndices]
            self.expMat = sparse_csr_matrix(data[1:, sampledCellIndices])
        else:
            self.geneList, self.sampleList, self.expMat = (
                data[1:, 0],
                data[0, 1:],
                sparse.csr_matrix(data[1:, 1:].astype("float")),
            )
        self.geneList = list(self.geneList)
        self.sampleList = list(self.sampleList)
        print("Genelist: " + str(self.geneList))
        del data
        self.expMat.eliminate_zeros()
        self.binMat = preprocessing.binarize(
            self.expMat, threshold=binarizeThreshold, copy=False
        )
        self.maxNodes = maxNodes
        self.maxSamples = maxSamples
        # self.binMat=self.expMat
        # super().__init__(self)
        # print(self.binMat.toarray())
        # populate binMat to a predefined size
        self.binMat.resize((min(len(self.geneList), 20000), len(self.sampleList)))
        self.pathwayGraphs = {}
    
    def __sampleCells(self, data, number_cells):
        """Sample a representative population of cells for rule inference - reduce memory requirements"""
        combined = np.apply_along_axis(lambda x: ''.join(str(x)), axis=1, arr=data)
        combined_weights = np.unique(combined, return_counts=True)[1]/len(combined)
        sampled_cells = np.random.choice(range(1, len(data[0, 1:])), size = number_cells, p = combined_weights)
        return(sampled_cells)

    def __addSubpop(self, subpopFile, sep):
        """Add subpopulation information to object"""
        # if(isinstance(subpop,list)):
        #    self.subpopInfo=zip(self.sampleList, subpopList)
        # else:
        #    if(isinstance(subpop, dict)):
        #        self.subpopInfo=subpops
        subpop_dict = {}
        subpopFile = open(subpopFile, "r")
        lines = subpopFile.readlines()
        for line in lines:
            cell = line.split(sep)
            subpop_dict[cell[0]] = cell[1].rstrip(
                "\n"
            )  # dict where key=cell and value=subpop
        self.subpopInfo = subpop_dict

    def __retrieveSubpop(self, cell):
        """Given a cell, retrieve the subpopulation information"""
        if hasattr(self, "subpopInfo"):
            return self.subpopInfo[cell]
        else:
            print("Please add subpopulation information first.")
            return None

    def __filterData(self, threshold):
        """Filter data based on CV cutoffs"""
        self.cvGenes = []
        if threshold is not None:
            for i in range(0, self.expMat.get_shape()[0]):
                rowData = list(self.expMat.getrow(i).todense())
                if np.std(rowData) / np.mean(rowData) >= threshold:
                    self.cvGenes.append(self.geneList[i])
        else:
            self.cvGenes = copy.deepcopy(self.geneList)

    def __sif_to_digraph(pathwayGenes, sif):
        node1, node2, edges = set(sif[0:, 0]), set(sif[0:, 2]), sif[:, [0, 2]]
        nodes = node1.union(node2)
        del node1, node2
        nodes = list(nodes)
        G = nx.DiGraph()
        G.add_nodes_from(nodes)
        for edge in sif:
            node1 = str(edge[0])
            node2 = str(edge[2])
            if node1 == node2:
                # don't add self edges
                continue
            elif (node1 not in pathwayGenes) or (node2 not in pathwayGenes):
                # don't add nodes that are not in the set of filtered genes
                continue
            elif node1 == node2:
                # don't add self-loops
                continue
            else:
                if int(edge[1]) == 1:
                    act = "a"
                elif int(edge[1]) == -1:
                    act = "i"
                else:
                    act = "u"
                G.add_edge(
                    str(edge[0]),
                    str(edge[2]),
                    signal=int(edge[1]),
                    activity=act,
                    interaction=act,
                )
        # graph post-processing
        # remove singletons/isolates
        G.remove_nodes_from(list(nx.isolates(G)))
        # To do: remove complexes, remove dependences of a node on complexes that include that node (which is a form of self-loop)
        return G

    def __makeSIF(pathway, keggObject):
        print("Pathway: ")
        print(pathway)
        activationRelations = [
            "activation",
            "binding/association",
            "phosphorylation",
            "indirect effect",
            "dissociation",
        ]  # Change
        inhibitionRelations = [
            "inhibition",
            "dephosphorylation",
            "dissociation",
            "ubiquitination",
        ]  # Change
        res = keggObject.parse_kgml_pathway(pathway)  # Change
        sif = []
        for rel in res["relations"]:
            # types can be PPrel (protein-protein interaction only)

            if rel["link"] != "PPrel":
                continue
            elif any(x in rel["name"] for x in activationRelations):
                Id1 = rel["entry1"]
                Id2 = rel["entry2"]
                type1 = res["entries"][[x["id"] for x in res["entries"]].index(Id1)][
                    "type"
                ]
                type2 = res["entries"][[x["id"] for x in res["entries"]].index(Id2)][
                    "type"
                ]
                if type1 != "gene" or type2 != "gene":
                    continue
                name1 = xstr(
                    res["entries"][[x["id"] for x in res["entries"]].index(Id1)][
                        "gene_names"
                    ]
                ).split(",")[0]
                name2 = xstr(
                    res["entries"][[x["id"] for x in res["entries"]].index(Id2)][
                        "gene_names"
                    ]
                ).split(",")[0]
                sif.append([name1, 1, name2])
            elif any(x in rel["name"] for x in inhibitionRelations):
                Id1 = rel["entry1"]
                Id2 = rel["entry2"]
                type1 = res["entries"][[x["id"] for x in res["entries"]].index(Id1)][
                    "type"
                ]
                type2 = res["entries"][[x["id"] for x in res["entries"]].index(Id2)][
                    "type"
                ]
                if type1 != "gene" or type2 != "gene":
                    continue
                name1 = xstr(
                    res["entries"][[x["id"] for x in res["entries"]].index(Id1)][
                        "gene_names"
                    ]
                ).split(",")[0]
                name2 = xstr(
                    res["entries"][[x["id"] for x in res["entries"]].index(Id2)][
                        "gene_names"
                    ]
                ).split(",")[0]
                sif.append([name1, -1, name2])
            else:
                pass
        sif = np.array([np.array(sif1) for sif1 in sif])
        return sif

    def __find_pathways(
        self, organism="hsa", minOverlap=20, writeGraphml=True, predefinedList=[]
    ):
        if hasattr(self, "cvGenes"):
            pathwayGenes = set(self.cvGenes)
        elif not hasattr(self, "cvGenes"):
            print("You have not filtered genes by any criterion.")
            pathwayGenes = set(self.geneList)
        k = KEGG()
        k.organism = organism
        if len(predefinedList) == 0:
            for x in list(k.pathwayIds):
                x = x.replace("path:", "")
                x = x.replace("hsa", "ko")
                try:
                    sif = makeSIF(pathway=x, keggObject=k)
                    if len(sif) < 1:
                        next
                    else:
                        nodes = set(sif[0:, 0]).union(set(sif[0:, 2]))
                        test = len(nodes.intersection(pathwayGenes))
                        if test >= minOverlap:
                            print(
                                "Pathway: ", x, " Overlap: ", test, " Edges: ", len(sif)
                            )
                            G = sif_to_digraph(sif=sif, pathwayGenes=pathwayGenes)
                            self.pathwayGraphs[x] = G
                            if writeGraphml:
                                nx.write_graphml(
                                    G,
                                    x + "_processed.graphml",
                                    infer_numeric_types=True,
                                )
                except:
                    next
        else:
            for x in list(predefinedList):
                # x=x.replace("path:","")
                # x=x.replace("hsa","ko")
                try:
                    sif = makeSIF(pathway=x, keggObject=k)
                    if len(sif) < 1:
                        next
                    else:
                        nodes = set(sif[0:, 0]).union(set(sif[0:, 2]))
                        test = len(nodes.intersection(pathwayGenes))
                        if test >= minOverlap:
                            print(
                                "Pathway: ", x, " Overlap: ", test, " Edges: ", len(sif)
                            )
                            G = sif_to_digraph(sif=sif, pathwayGenes=pathwayGenes)
                            self.pathwayGraphs[x] = G
                            if writeGraphml:
                                nx.write_graphml(
                                    G,
                                    x + "_processed.graphml",
                                    infer_numeric_types=True,
                                )
                except:
                    next

    def __add_pathways(
        self, pathwayList=[], minOverlap=20, writeGraphml=True, removeSelfEdges=False
    ):
        """Add a list of pathways in graphml format to the singleCell object"""
        if hasattr(self, "cvGenes"):
            pathwayGenes = set(self.cvGenes)
        elif not hasattr(self, "cvGenes"):
            print("You have not filtered genes by any criterion.")
            pathwayGenes = set(self.geneList)

        if isinstance(pathwayList, list):
            for (
                x
            ) in (
                pathwayList
            ):  # list(glob.glob("*.graphml")):  # list(glob.glob('*[0-9].graphml')):
                G = nx.read_graphml(x)
                nodes = set(G.nodes())
                test = len(nodes.intersection(pathwayGenes))

                if test >= minOverlap:
                    print(
                        "Pathway: ", x, " Overlap: ", test, " Edges: ", len(G.edges())
                    )
                    nodes = list(G.nodes())
                    if removeSelfEdges:
                        G.remove_edges_from(nx.selfloop_edges(G))  # remove self loops
                    # remove genes not in dataset
                    for pg in list(G.nodes()):
                        if pg not in pathwayGenes:
                            G.remove_node(pg)
                    # graph post-processing
                    # remove singletons/isolates
                    G.remove_nodes_from(list(nx.isolates(G)))
                    # To do: remove complexes, remove dependences of a node on complexes that include that node (which is a form of self-loop)
                    self.pathwayGraphs[x] = G
                    print(
                        "Edges after processing:",
                        len(G.edges()),
                        " Overlap: ",
                        len(set(G.nodes()).intersection(pathwayGenes)),
                    )
                    if writeGraphml:
                        nx.write_graphml(
                            G, x + "_processed.graphml", infer_numeric_types=True
                        )
        else:
            if isinstance(pathwayList, dict):
                for x, G in pathwayList.items():
                    nodes = set(G.nodes())
                    test = len(nodes.intersection(pathwayGenes))

                    if test >= minOverlap:
                        print(
                            "Pathway: ",
                            x,
                            " Overlap: ",
                            test,
                            " Edges: ",
                            len(G.edges()),
                        )
                        nodes = list(G.nodes())
                        if removeSelfEdges:
                            G.remove_edges_from(
                                nx.selfloop_edges(G)
                            )  # remove self loops
                        # remove genes not in dataset
                        for pg in list(G.nodes()):
                            if pg not in pathwayGenes:
                                G.remove_node(pg)
                        # graph post-processing
                        # remove singletons/isolates
                        G.remove_nodes_from(list(nx.isolates(G)))
                        # To do: remove complexes, remove dependences of a node on complexes that include that node (which is a form of self-loop)
                        self.pathwayGraphs[x] = G
                        print(
                            "Edges after processing:",
                            len(G.edges()),
                            " Overlap: ",
                            len(set(G.nodes()).intersection(pathwayGenes)),
                        )
                        if writeGraphml:
                            nx.write_graphml(
                                G, x + "_processed.graphml", infer_numeric_types=True
                            )

    def __makemetaNetwork(self):
        if not hasattr(self, "pathwayGraphs"):
            print("Run find_pathways before trying to construct a meta network")
        else:
            graphs = list(self.pathwayGraphs.values())
            self.metaNetwork = nx.compose_all(graphs)
            # for item in nx.strongly_connected_components(self.metaNetwork):
            #     print(item)
            largest = max(nx.strongly_connected_components(self.metaNetwork), key=len)
            # largest=max(nx.kosaraju_strongly_connected_components(self.metaNetwork), key=len) #same result in trials
            self.metaNetwork.remove_nodes_from(
                [n for n in self.metaNetwork if n not in list(largest)]
            )
            print(
                "The meta network has ", self.metaNetwork.number_of_nodes(), " nodes."
            )

    def __genInitValueList(self, graph):

        newInitValueList = []
        numberCells = len(self.sampleList)
        nodes = list(graph)
        numberNodes = len(nodes)

        geneIndex = []
        for g in self.geneList:
            if g in nodes:
                i = self.geneList.index(g)
                geneIndex.append(i)

        # print(geneIndex)

        for j in range(0, numberCells):
            temp = self.binMat[geneIndex, j].todense().tolist()
            temp = [y for x in temp for y in x]
            newInitValueList.append(temp)

        # print(len(newInitValueList))

        # return newInitValueList
        self.initValueList = newInitValueList

    def __setupEmptyKOKI(self):
        """make empty list representing no knockouts or knockins"""
        self.knockoutLists = [0] * len(self.nodePositions)
        self.knockinLists = [0] * len(self.nodePositions)

    def __scoreNodes(self, graph):
        """Wrapper function - performs rule determination and node scoring in preparation for pathway analysis"""

        net = nx.read_graphml(graph)
        # print(self.geneList)
        # print(list(net))
        netGenes = [
            self.geneList.index(gene) for gene in list(net) if gene in self.geneList
        ]
        print(netGenes)
        updateBooler = ctypes.cdll.LoadLibrary("./simulator.so")
        scSyncBoolC = updateBooler.scSyncBool
        importanceScore = updateBooler.importanceScore
        self._singleCell__inherit(
            net,
            removeSelfEdges=False,
            restrictIncomingEdges=True,
            maxIncomingEdges=3,
            groundTruth=False,
        )
        self._ruleMaker__updateCpointers()
        self.__setupEmptyKOKI()

        # Genetic algorithm
        population, logbook = self._ruleMaker__eaMuPlusLambdaAdaptive(
            scSyncBoolC, graph
        )
        out1, out2, model = self._ruleMaker__findPopBest(population)
        with open(graph + "_rules_GA.txt", "w") as text_file:
            text_file.write(model.writeModel(out2, model))
        pickle.dump(out2, open(graph + "_out2.pickle", "wb"))

        # Local search
        outputs = [
            self._ruleMaker__checkNodePossibilities(
                node, out2, model.knockoutLists, model.knockinLists, scSyncBoolC
            )
            for node in range(0, len(model.nodePositions))
        ]
        equivs = []
        individual = []
        devs = []
        localErrors = []
        for output in outputs:
            individual.extend(output[0])
            equivs.append(output[1])
            devs.append(output[2])
            localErrors.append(output[3])
        pickle.dump(equivs, open(graph + "_equivs1.pickle", "wb"))
        pickle.dump(localErrors, open(graph + "_localErrors1.pickle", "wb"))
        bruteout2 = []
        for i in range(len(equivs)):
            bruteout2.extend(equivs[i][randint(0, len(equivs[i]) - 1)])
        with open(graph + "_rules_LS.txt", "w") as text_file:
            text_file.write(model.writeModel(bruteout2, model))
        pickle.dump(bruteout2, open(graph + "_bruteout2.pickle", "wb"))

        # Importance score calculation
        importanceScores = self._ruleMaker__calcImportance(
            "", self, importanceScore, graph
        )
        print(importanceScores)

    def __scorePathway(self, RAs, pathImportances, pathname):
        """calculate z score for a given pathway"""
        CVdict = {}
        impactScore = {}
        score = 0
        allNodes = list(RAs.keys())
        for node in pathImportances:
            print("***Node: " + str(node) + "***")
            print("RA: " + str(RAs[node]))
            print("IS: " + str(pathImportances[node]))
            nodeData = (
                self.binMat[self.geneList.index(node), range(0, len(self.sampleList))]
                .todense()
                .tolist()
            )
            # score+=abs(float(RAs[node]))*math.log(float(pathImportances[node]),2) #*CVdict[node]
            # score+=abs(float(RAs[node]))*float(pathImportances[node])
            CVdict[node] = np.std(nodeData)
            print(CVdict[node])
            impactScore[node] = (
                abs(float(RAs[node])) * float(pathImportances[node]) * CVdict[node]
            )
            # print(impactScore[node])
            score += impactScore[node]
            # score+=(abs(float(RAs[node]))*float(pathImportances[node]))/nodeData_std
        # print(score)
        print(allNodes)
        print(
            "Relative abundance mean difference: "
            + str(np.mean([abs(RAs[value]) for value in allNodes]))
        )
        randomScores = []
        for i in range(1000):
            tempscore = 0
            for node in pathImportances:
                t1 = str(allNodes[randint(0, len(allNodes) - 1)])
                # tempscore+=abs(RAs[t1]*math.log(float(pathImportances[node]),2)) #*CVdict[node]
                # tempscore+=abs(RAs[t1]*float(pathImportances[node]))
                tempscore += abs(RAs[t1] * float(pathImportances[node])) * CVdict[node]
                # tempscore+=(abs(RAs[t1]*float(pathImportances[node])))/CVdict[node]
            randomScores.append(tempscore)
        meaner = np.mean(randomScores)
        stdev = np.std(randomScores)
        zscore = (score - meaner) / stdev

        # make histogram of the scores
        # An "interface" to matplotlib.axes.Axes.hist() method
        n, bins, patches = plt.hist(
            x=randomScores, bins="auto", color="#0504aa", alpha=0.7, rwidth=0.85
        )
        plt.grid(axis="y", alpha=0.75)
        plt.xlabel("Impact Score")
        plt.ylabel("Frequency")
        maxfreq = n.max()
        # Set a clean upper y-axis limit.
        plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
        plt.title(self.__getPathwayName(str(pathname)))
        axes = plt.gca()
        y_min, y_max = axes.get_ylim()
        plt.text(
            meaner,
            y_max - 15,
            "Mean="
            + str(round(meaner, 2))
            + "\n Std.dev="
            + str(round(stdev, 2))
            + "\n Pathway score="
            + str(round(score, 2)),
            bbox=dict(facecolor="red", alpha=0.75),
        )
        plt.savefig(
            pathname[:8] + "_impact_scores_hist.png", format="png", bbox_inches="tight"
        )
        plt.close()
        return zscore, impactScore, CVdict

    def run_pathway_analysis(
        self,
        contrast,
        conditions,
        conditionsSep="\t",
        contrastSep="\t",
        pathwayList=None,
    ):

        """
        Uses the output files from rule determination to perform pathway analysis for either specified networks or all the networks initially used for rule inference.

        Parameters:
        contrast: path to a text file where each line contains the two conditions (corresponding to column labels in the conditions file) are to be compared during pathway analysis.

        conditions: path to the file containing cell-wise condition labels, ie, metadata for each cell in the training dataset. The columns are condition   variables and the rows are cells. The first column must be cell names that correspond to the columns in the training data file. The column names must contain the variables specified in the contrast file

        conditions_separator: Separator for the conditions file. Must be one of 'comma', 'space' or 'tab' (spell out words, escape characters will not work)."

        contrast_separator: Separator for the contrast file. Must be one of 'comma', 'space' or 'tab' (spell out words, escape characters will not work)."

        pathwayList: Paths to GRAPHML files that should be used for scBONITA analysis. Usually networks from non-KEGG sources, saved in GRAPHML format
        The default value is to use sll the networks initially used for rule inference.
        """
        if pathwayList is None:
            pathwayList = list(self.pathwayGraphs.keys())
        self.binMat2 = self.binMat.A
        self.expMat2 = self.expMat.A
        pvalDict = {}
        overallUpreg = {}
        # open the set of differences to be considered
        contrastList = []
        contrast = pd.read_csv(contrast, sep=contrastSep, header=None)
        # for row in csv.reader(open(contrast), delimiter=contrastSep):
        #    contrastList.append(list(row))
        for index, row in contrast.iterrows():
            contrastList.append([row[0], row[1]])
        print(contrastList)
        conditions = pd.read_csv(conditions, sep=conditionsSep, index_col=0)
        print(conditions)
        for contrasts in contrastList:
            print(contrasts)
            # condition1 = conditions.loc[:, contrasts[0]] #[[str(contrasts[0])]]
            # condition2 = conditions.loc[:, contrasts[1]]#[[str(contrasts[1])]]
            cells_condition1 = list(
                conditions.index[conditions.loc[:, contrasts[0]] == 1]
            )
            cells_condition2 = list(
                conditions.index[conditions.loc[:, contrasts[1]] == 1]
            )
            index_condition1 = [
                self.sampleList.index(i)
                for i in set(cells_condition1).intersection(set(self.sampleList))
            ]
            index_condition2 = [
                self.sampleList.index(i)
                for i in set(cells_condition2).intersection(set(self.sampleList))
            ]
            # make RA - in the case of single cell experiments, find the proportion of cells in which the gene is expressed
            # for pathname in list(self.pathwayGraphs.keys()):
            for pathname in pathwayList:
                ISfile = glob.glob(pathname + "*_importanceScores.csv")
                paNetTemp = glob.glob(pathname + "*_processed.graphml")
                if len(ISfile) > 0 and len(paNetTemp) > 0:
                    print(ISfile[0])
                    paNetTemp = nx.read_graphml(paNetTemp[0])
                    # get nodeScores
                    nodeScoresDF = pd.read_csv(ISfile[0])
                    nodeScoresDF.index = list(nodeScoresDF.Node)
                    if "Importance Score" in nodeScoresDF.columns:
                        # add impact score as attribute to graph
                        nodeScores = dict(
                            zip(
                                list(nodeScoresDF.index),
                                list(nodeScoresDF.loc[:, "Importance Score"]),
                            )
                        )
                        nx.set_node_attributes(
                            paNetTemp, values=nodeScores, name="importanceScore"
                        )
                        # add obsERS as attribute to graph
                        obsERS = dict(
                            zip(
                                list(nodeScoresDF.index),
                                list(nodeScoresDF.loc[:, "ObsERS"]),
                            )
                        )
                        nx.set_node_attributes(
                            paNetTemp, values=obsERS, name="Observed ERS"
                        )
                        # add maxERS as attribute to graph
                        maxERS = dict(
                            zip(
                                list(nodeScoresDF.index),
                                list(nodeScoresDF.loc[:, "MaxERS"]),
                            )
                        )
                        nx.set_node_attributes(paNetTemp, values=obsERS, name="Max ERS")
                        # write out IS graph with additions
                        # Calculate RA
                        RA = {}
                        upreg_condition1 = {}
                        expression_condition1 = {}
                        expression_condition2 = {}
                        nodeScoresDF[str(contrasts[0])] = np.nan
                        nodeScoresDF[str(contrasts[1])] = np.nan
                        for node in list(nodeScoresDF.index):
                            node_index = self.geneList.index(node)
                            expression_condition1[node] = np.mean(
                                self.expMat2[node_index, index_condition1].tolist()
                            )
                            expression_condition2[node] = np.mean(
                                self.expMat2[node_index, index_condition2].tolist()
                            )
                            RA[node] = (
                                expression_condition1[node]
                                - expression_condition2[node]
                            )
                            if (
                                expression_condition1[node]
                                > expression_condition2[node]
                            ):
                                upreg_condition1[node] = True
                            else:
                                upreg_condition1[node] = False
                        nodeScoresDF[str(contrasts[0])] = nodeScoresDF["Node"].map(
                            expression_condition1
                        )
                        nodeScoresDF[str(contrasts[1])] = nodeScoresDF["Node"].map(
                            expression_condition2
                        )
                        nodeScoresDF[
                            str("Upregulated_in_" + str(contrasts[0]))
                        ] = nodeScoresDF["Node"].map(upreg_condition1)
                        # add RA as attribute to graph
                        nx.set_node_attributes(
                            paNetTemp, values=RA, name="relativeAbundance"
                        )
                        nx.set_node_attributes(
                            paNetTemp,
                            values=expression_condition1,
                            name=str(contrasts[0]),
                        )
                        nx.set_node_attributes(
                            paNetTemp,
                            values=expression_condition2,
                            name=str(contrasts[1]),
                        )
                        nx.set_node_attributes(
                            paNetTemp,
                            values=upreg_condition1,
                            name=str("Upregulated_in_" + str(contrasts[0])),
                        )
                        z_scores = []
                        modRA = {}
                        for node in list(nodeScoresDF.index):
                            modRA[node] = abs(RA[node])
                        # iterate over comparisons for each pathway and calculate z score
                        zscore, impactScore, CVdict = self.__scorePathway(
                            modRA, nodeScores, pathname[:-8]
                        )
                        z_scores.append(zscore)
                        pvals = scipy.stats.norm.sf(z_scores)  # calculate p value
                        # print(pvals)
                        pvalDict[str(pathname)] = [
                            pvals,
                            # str(CVdict),
                            # str(zscore),
                        ]
                        print(pvalDict[str(pathname)])
                        # add impact score as attribute to graph
                        # nx.set_node_attributes(paNetTemp, values=impactScore, name='impactScore')

                        # write out graph with additions if pval < 0.05
                        if pvalDict[str(pathname)][0] < 1:
                            nx.write_graphml_lxml(
                                paNetTemp,
                                pathname + "_IS_" + "_vs_".join(contrasts) + ".graphml",
                            )
                        nodeScoresDF.to_csv(
                            pathname
                            + ".graphml_processed.graphml_importanceScores.csv",
                            index=False,
                        )
                        if (
                            sum(
                                nodeScoresDF[str("Upregulated_in_" + str(contrasts[0]))]
                            )
                            > len(list(nodeScoresDF.index)) / 2
                        ):
                            overallUpreg[str(pathname)] = str("True")
                        else:
                            overallUpreg[str(pathname)] = str("False")
                else:
                    print("No output for " + pathname)
                    # pvalDict[str(pathname)]=["None"]
                    # overallUpreg[str(pathname)] = ["None"]

            self.pvalDict = pvalDict
            # fh=open("pvalues_"+"_vs_".join(contrasts)+".csv", "w+")
            # fh.write(','.join(["Pathway ID", "Pathway Name", "P value", "Contrast", str("Upregulated_in_"+str(condition1)+"\n")]))
            for key, value in pvalDict.items():
                key = str(key)
                pvalDict[key] = [
                    key,
                    str(self.__getPathwayName(key)),
                    str(value[0][0]),
                    "_vs_".join(contrasts),
                    str(overallUpreg[key]),
                    # str(value[1]),
                    # str(value[2])
                ]
            pvalDF = pd.DataFrame.from_dict(pvalDict, orient="index")
            pvalDF.columns = [
                "Pathway ID",
                "Pathway Name",
                "P value",
                "Contrast",
                str("Upregulated_in_" + str(contrasts[0])),
            ]
            pvalDF.loc[:, "P value"] = [
                float(temp) for temp in pvalDF.loc[:, "P value"]
            ]
            print(pvalDF.loc[:, "P value"])
            pvalDF.loc[:, "Adjusted P value"] = list(
                multipletests(list(pvalDF.loc[:, "P value"]), method="bonferroni")[1]
            )
            pvalDF.to_csv("pvalues_" + "_vs_".join(contrasts) + ".csv", index=False)
        del self.binMat2

    def __inherit(
        self,
        graph,
        removeSelfEdges=False,
        restrictIncomingEdges=True,
        maxIncomingEdges=3,
        groundTruth=False,
        graphName="",
    ):
        super().__init__(
            graph,
            removeSelfEdges,
            restrictIncomingEdges,
            maxIncomingEdges,
            groundTruth,
            graphName,
        )

    def __getPathwayName(self, hsaURL):
        fileReg = re.compile("NAME\s+(\w+.*)")
        pathwayFile = requests.get("http://rest.kegg.jp/get/" + hsaURL, stream=True)
        for line in pathwayFile.iter_lines():
            line = line.decode("utf-8")
            result = fileReg.match(line)
            if result:
                return result.group(1)
        return hsaURL

    def processERS_minimalRule(self, equivsName):
        """Create an individual from the ERS generated by the local search, for importance score calculation"""
        ersFile = open(str(equivsName), "rb")
        ers = pickle.load(ersFile)
        ersFile.close()
        individual = []
        for i in range(len(ers)):
            ruleSet = ers[i]
            numOrRules = [sum(ruleSet[j]) for j in range(len(ruleSet))]
            deciderVar = max(numOrRules)  # maximal or rules
            maxOrRules = ruleSet[
                numOrRules.index(deciderVar)
            ]  # rules with maximum or terms
            maxUpstreamNodes = 0
            minimalRule = []
            for orRule in [maxOrRules]:
                if sum(orRule) > 0:
                    numUpstreamNodes = [
                        self.andNodeList[i][orTerm]
                        for orTerm in range(len(orRule))
                        if orRule[orTerm] == 1
                    ]
                else:
                    minimalRule = orRule
                    continue
                numUpstreamNodes = [len(element) for element in numUpstreamNodes]
                numUpstreamNodes = sum(numUpstreamNodes)
                if numUpstreamNodes > maxUpstreamNodes:
                    maxUpstreamNodes = numUpstreamNodes
                    minimalRule = orRule
                else:
                    maxUpstreamNodes = maxUpstreamNodes
                    minimalRule = minimalRule
            individual.extend(minimalRule)
        return individual

    def processERS(self, equivsName):
        """Create an individual from the ERS generated by the local search, for importance score calculation"""
        ersFile = open(str(equivsName), "rb")
        ers = pickle.load(ersFile)
        ersFile.close()
        individual = []
        for i in range(len(ers)):
            ruleSet = ers[i]
            rand = random.randint(0, len(ruleSet) - 1)
            randomRule = ruleSet[rand]
            individual.extend(randomRule)
        return individual

    def __updateBool(
        self,
        currentNode,
        oldValue,
        nodeIndividual,
        andNodes,
        andNodeInvertList,
        nodeStart,
        nodeEnd,
    ):
        indindex = nodeStart
        orset = [np.nan for i in range(self.maxNodes)]
        counter = 0
        while indindex < nodeEnd:
            andindex = indindex - nodeStart
            if nodeIndividual[indindex]:
                newval = (
                    oldValue[andNodes[andindex][0]] != andNodeInvertList[andindex][0]
                )
                addnode = 0
                while addnode < len(andNodes[andindex]):
                    if andNodes[andindex][addnode] > (-1):
                        newval = newval and (
                            oldValue[andNodes[andindex][addnode]]
                            != andNodeInvertList[andindex][addnode]
                        )
                    addnode = addnode + 1
                orset[counter] = newval
                counter = counter + 1
            indindex = indindex + 1
        newval = orset[0]
        q = 1
        while q < counter:
            newval = newval or orset[q]
            q = q + 1
        return newval

    def __getAttractCycle2(self, valsT):
        nodeNum = len(self.nodeList)
        STEP = 10
        res = [-1, -1]
        i = STEP - 1
        while i >= 0:
            flag = 0
            j = i - 1
            while j >= 0:
                flag = 1
                k = 0
                while k < nodeNum:
                    if valsT[i][k] != valsT[j][k]:
                        flag = 0
                    k = k + 1
                if flag == 1:
                    res[0] = j
                    res[1] = i
                    return res
                j = j - 1
            i = i - 1
        return res

    def __getAttract(self, vals, resSubmit):
        res = self.__getAttractCycle2(vals)
        resSubmit[0] = res[0]
        resSubmit[1] = res[1]
        return resSubmit

    def __cluster(
        self,
        simData,
        resSubmit,
        sampleIndex,
        individual,
        indLen,
        nodeNum,
        andLenList,
        individualParse,
        andNodes,
        andNodeInvertList,
        simSteps,
        knockouts,
        knockins,
        binMat,
        nodePositions,
    ):
        gc.collect()
        oldValue = [np.nan for i in range(0, nodeNum)]
        newValue = [np.nan for i in range(0, nodeNum)]
        newValue = [binMat[nodePositions[i]][sampleIndex] for i in range(nodeNum)]
        simData[0] = [newValue[i] for i in range(nodeNum)]
        step = 1
        while step < simSteps:
            oldValue = [newValue[i] for i in range(nodeNum)]
            i = 0
            while i < nodeNum:
                if knockouts[i] == 1:
                    temp = 0
                    newValue[i] = temp
                    simData[step][i] = temp
                elif knockins[i] == 1:
                    temp = 1
                    newValue[i] = temp
                    nodePos = nodePositions[i]
                    simData[step][i] = temp
                elif andLenList[i] == 1:
                    temp = oldValue[andNodes[i][0][0]] != andNodeInvertList[i][0][0]
                    newValue[i] = temp
                    nodePos = nodePositions[i]
                    simData[step][i] = temp
                elif andLenList[i] == 0:
                    temp = oldValue[i]
                    newValue[i] = temp
                    nodePos = nodePositions[i]
                    simData[step][i] = temp
                elif i == (nodeNum - 1):
                    nodeEnd = indLen
                else:
                    nodeEnd = individualParse[i + 1]
                    nodeStart = individualParse[i]
                    temp = self.__updateBool(
                        i,
                        oldValue,
                        individual,
                        andNodes[i],
                        andNodeInvertList[i],
                        nodeStart,
                        nodeEnd,
                    )
                    newValue[i] = temp
                    nodePos = nodePositions[i]
                    simData[step][i] = temp
                i = i + 1
            step = step + 1
        gc.collect()
        resSubmit = self.__getAttract(simData, resSubmit)
        gc.collect()
        return resSubmit

    def assignAttractors(
        self,
        pathwayFiles=[],
        useMinimalRuleSet=True,
        simSteps=10,
        removeZeroAttractors=True,
    ):

        """Uses a dataframe of cells and the attractors to which they are assigned (usually the output of assignAttractors) to generate (a) a barplot showing the frequency of cells assigned to (selected) attractors, optionally faceted by user-specified variables and (b) a 2D scatterplot, such as a UMAP or tSNE plot, showing cells colored by the attractors to which they are assigned.

        Parameters
        ----------

        plottingData: pandas DataFrame
            Usually a 2-column dataframe where rows = cells (indices must contain the cells in the training dataset) and columns = UMAP or tSNE (or similar) dimensions
        distanceDF: pandas DataFrame
            Output of assignAttractors, usually. A pandas DataFrame where rows = cells (indices must contain the cells in the training dataset) and at least one column (specified by attractorsAttribute parameter, see below) contains the attractors to which the cells have been assigned
        attractorsAttribute: str
            The column in distanceDF which contains the attractors to which the cells have been assigned. Default is 'decider', as this is the name in the output of assignAttractors
        allAttractors: bool
            True if all assigned attractors are to be included in the frequency and UMAP plots. Should usually be 'True' for making preliminary analysis plots, and 'False' while making publication-quality plots, for ease of visualization and to leave out infrequent attractors.
        numberOfAttractorsToShow: int
            The number of attractors to be shown in the analysis plots. Used when allAttractors is False. Only the top (numberOfAttractorsToShow) most frequent attractors will be shown. The others will be collapsed into a single category "Others".
        cmap: str
            The matplotlib color palette to be used for the analysis plots. Consider using 'colorblind' if the number of categories/attractors is small and a continuous palette such as 'Blues' or 'Greens' if the number is large. See matplotlib documentation for more options on color palettes.
        makeFrequencyPlots: bool
            Should frequency plots be made?
        frequencyGrouping: list
            List of variables in the plottingData which should be used for faceting the frequency plots. Eg. sample, batch, disease/control, etc.
        freqplotsFile: str
            Path to output PDF file for frequency plots. Default plots to the standard output.
        makeUMAP: bool
            Should UMAP plots be made?
        umapFile: str
            Path to output PDF file for UMAP plots. Default plots to the standard output.

        Returns
        -------
            None
        """

        # get the networks used for rule inference - you can change this parameter when you call the function
        if len(pathwayFiles) == 0:
            pathwayFiles = glob.glob("*_processed.graphml")
        else:
            pathwayFiles = pathwayFiles
        # Reconstruct the singleCell object
        # objectFile = glob.glob("*.binscTest.pickle")
        # scObject = pickle.load(open(objectFile[0], "rb"))
        featureTable = {}  # record of assignment of cells to attractors
        enumerateDict = {}  # records attractor counts by network
        # iterate over all networks used for rule inference
        for pathway in pathwayFiles:
            # if path.exists(str(pathway + "_equivs1.pickle")) and not path.exists(
            #        str(pathway + "_attractorDistance_mgp.csv")
            # ):  #if rules were successfully inferred and assignAttractors has not yet been successfully executed
            if path.exists(str(pathway + "_equivs1.pickle")):
                print("#############")
                # Read in graph
                print(pathway)
                net = nx.read_graphml(pathway)
                # Create ruleMaker object
                self.__inherit(
                    net,
                    removeSelfEdges=False,
                    restrictIncomingEdges=True,
                    maxIncomingEdges=3,
                    groundTruth=False,
                )
                self.__setupEmptyKOKI()
                self._ruleMaker__updateCpointers()
                KOlist = []  # knocked-out genes
                KIlist = []  # knocked-in genes
                attractorList = []
                # individual = processERS(
                #    str(pathway + "_equivs1.pickle"))  #get a rule set #
                if useMinimalRuleSet:
                    individual = self.processERS_minimalRule(
                        str(pathway + "_equivs1.pickle")
                    )  # get a rule set
                else:
                    individual = self.processERS(
                        str(pathway + "_equivs1.pickle")
                    )  # get a rule set
                with open(pathway + "_minimal_rules.txt", "w") as text_file:
                    text_file.write(self.writeModel(individual, self))
                text_file.close()
                knockins = np.zeros(
                    len(self.nodeList), dtype=np.intc, order="C"
                )  # set up empty ko and ki genes
                knockouts = np.zeros(len(self.nodeList), dtype=np.intc, order="C")
                for knocker in KOlist:
                    knockouts[knocker] = 1
                for knocker in KIlist:
                    knockins[knocker] = 1
                # put objects in correct format for passing to C
                nodeIndividual = np.array(individual, dtype=np.intc, order="C")
                indLen = len(nodeIndividual)
                andNodes = self.andNodeList
                nodeNum = len(self.nodeList)
                andNodeInvert = self.andNodeInvertList
                individualParse = np.array(
                    self.individualParse, dtype=np.intc, order="C"
                )
                andLenList = np.array(self.andLenList, dtype=np.intc, order="C")
                nodePositions1 = self.nodePositions
                nodePositionsC = np.array(nodePositions1, dtype=np.intc, order="C")
                binMatC3 = np.array(
                    copy.deepcopy(self.binMat.toarray(order="C")),
                    order="C",
                    dtype=np.intc,
                )
                for sampleIndex in range(0, len(self.sampleList)):
                    vals = np.full(
                        shape=(simSteps, nodeNum),
                        fill_value=0,
                        dtype=np.intc,
                        order="C",
                    )
                    res = np.full(
                        shape=2, fill_value=2, dtype=np.intc, order="C"
                    )  # initiate output array
                    res = self._singleCell__cluster(
                        vals,
                        res,
                        sampleIndex,
                        nodeIndividual,
                        indLen,
                        nodeNum,
                        andLenList,
                        individualParse,
                        andNodes,
                        andNodeInvert,
                        simSteps,
                        knockouts,
                        knockins,
                        binMatC3,
                        self.nodePositions,
                    )
                    attractor = []
                    for temp in range(res[0], res[1] + 1):
                        attractor.append(tuple(vals[temp, :]))
                    attractorList.extend(attractor)
                # Find unique attractor states
                attractorList = [list(x) for x in set(tuple(x) for x in attractorList)]
                print(
                    str(len(attractorList))
                    + " unique attractors were identified for "
                    + pathway
                )
                if removeZeroAttractors:
                    # remove attractors that are all 0
                    testIfZero = [sum(x) for x in set(tuple(x) for x in attractorList)]
                    attractorList = [
                        attractorList[i]
                        for i in range(len(attractorList))
                        if testIfZero[i] > 0
                    ]
                attractorList = tuple(attractorList)
                # Calculate new Distance between chosen attractors and cells
                newdistanceDict = {}
                for i in range(0, len(self.sampleList)):
                    cell = str(self.sampleList[i])
                    # get value of cell
                    cellInitValue = self.binMat[:, i]
                    cellInitValue = cellInitValue.toarray()
                    cellInitValue = cellInitValue.tolist()
                    cellInitValue = [
                        cellInitValue[temp] for temp in list(nodePositions1)
                    ]
                    # set up data structure to store Hamming distance to attractors
                    newdistanceDict[cell] = []
                    for attr in attractorList:
                        newdistanceDict[cell].append(
                            sum(
                                [abs(i1 - i2[0]) for i1, i2 in zip(attr, cellInitValue)]
                            )
                        )
                # Change table to a pandas dataframe for easy readability and plotting
                distanceDF = pd.DataFrame.from_dict(newdistanceDict, orient="index")
                distanceDF.columns = [str(attr) for attr in attractorList]
                enumerateDict[pathway] = str(len(attractorList))
                # distanceDF['decider'] = [
                #    np.inf for i in range(0, len(distanceDF.index))
                # ]
                decider = distanceDF.idxmin(axis=1)
                decider = [distanceDF.columns.get_loc(temp) for temp in decider]
                distanceDF.loc[:, "deciderDistance"] = distanceDF.min(axis=1)
                distanceDF.loc[:, "decider"] = decider
                print(
                    str(len(set(decider)))
                    + " of these unique attractors mapped to cells for "
                    + pathway
                )
                distanceDF.to_csv(str(pathway + "_attractorDistance.csv"), index=True)
                return distanceDF
            else:
                print("No output from rule inference for " + pathway)

    def makeAttractorHeatmaps(
        self,
        distanceDF,
        network,
        width,
        height,
        allAttractors=True,
        cmap="RdYlBu_r",
        numberOfAttractorsToShow=False,
        outputFileName=None,
    ):
        sns.set_context("paper", font_scale=0.5)
        numberClusters = len(distanceDF.decider.unique())
        attractorList = pd.Series(distanceDF.columns)[0 : (len(distanceDF.columns) - 2)]
        attractorList = attractorList.apply(literal_eval)
        attractorList = list(set([tuple(temp) for temp in attractorList]))
        attractorDF = pd.DataFrame(attractorList).T
        attractorDF.index = [self.geneList[temp] for temp in self.nodePositions]
        attractorDF = attractorDF.loc[:, attractorDF.columns.isin(distanceDF.decider)]
        plt.figure(figsize=(width, height))
        print("Attractor cell counts: ")
        for i in list(set(distanceDF.decider)):
            print(i, list(distanceDF.decider).count(i))

        attractorDF = attractorDF.T
        attractorDF = attractorDF[
            [i for i in attractorDF if len(set(attractorDF[i])) > 1]
        ]
        attractorDF = attractorDF.T
        attractorDF2 = attractorDF
        if isinstance(allAttractors, list):
            attractorDF2 = attractorDF.reindex(columns=allAttractors)
            attractorDF2 = attractorDF2.loc[attractorDF2.sum(axis=1) != 0, :]
            attractorDF2 = attractorDF2.loc[
                attractorDF2.sum(axis=1) != len(allAttractors), :
            ]
        else:
            if not allAttractors and numberOfAttractorsToShow:
                attractorDF2 = attractorDF
                allAttractors = distanceDF.decider.value_counts()[
                    :numberOfAttractorsToShow
                ].index.tolist()  # top n attractors where n = numberOfAttractorsToShow
                attractorDF2 = attractorDF.loc[:, allAttractors]
                attractorDF2 = attractorDF2.loc[attractorDF2.sum(axis=1) != 0, :]
                attractorDF2 = attractorDF2.loc[
                    attractorDF2.sum(axis=1) != numberOfAttractorsToShow, :
                ]
            else:
                if not allAttractors and not numberOfAttractorsToShow:
                    attractorDF2 = attractorDF
                    allAttractors = (
                        distanceDF.decider.value_counts().index.tolist()
                    )  # top n attractors where n = numberOfAttractorsToShow
                    attractorDF2 = attractorDF.loc[:, allAttractors]
                    attractorDF2 = attractorDF2.loc[attractorDF2.sum(axis=1) != 0, :]
                    attractorDF2 = attractorDF2.loc[
                        attractorDF2.sum(axis=1) != numberClusters, :
                    ]
        if attractorDF2.shape[1] > 1:
            ax = sns.clustermap(
                attractorDF2,
                cmap=cmap,
                linewidth=0.3,
                yticklabels=True,
                xticklabels=True,
                vmin=0,
                vmax=1,
                cbar_kws={"ticks": [0, 1]},
                figsize=(width, height),
            )
            ax.cax.set_position([0.15, 0.2, 0.03, 0.45])
            ax.cax.yaxis.set_label_position("left")
            ax.ax_col_dendrogram.set_visible(False)
            ax.ax_row_dendrogram.set_visible(False)
            plt.show()
            if outputFileName is not None:
                ax.savefig(outputFileName, bbox_inches="tight", pad_inches=0.01)
            else:
                ax.savefig(
                    str(network) + "_multiple_attr_heatmap.pdf",
                    bbox_inches="tight",
                    pad_inches=0.01,
                )
            plt.clf()
        else:
            ax = sns.heatmap(
                attractorDF2,
                vmin=0,
                vmax=1,
                cmap=cmap,
                xticklabels=False,
                cbar_kws={"ticks": [0, 1]},
            )
            ax.set(ylabel="Genes")
            if outputFileName is not None:
                ax.savefig(outputFileName, bbox_inches="tight", pad_inches=0.01)
            else:
                ax.figure.savefig(
                    str(network) + "_one_attr_heatmap_mgp.pdf",
                    bbox_inches="tight",
                    pad_inches=0.01,
                )
            plt.clf()
        return (attractorDF, ax)

    def makeAttractorAnalysisPlots(
        self,
        plottingData,
        distanceDF,
        attractorsAttribute="decider",
        allAttractors=False,
        numberOfAttractorsToShow=2,
        cmap="colorblind",
        makeFrequencyPlots=True,
        frequencyGrouping="",
        freqplotsFile="",
        makeUMAP=True,
        umapFile="",
    ):

        """Uses a dataframe of cells and the attractors to which they are assigned (usually the output of assignAttractors) to generate (a) a barplot showing the frequency of cells assigned to (selected) attractors, optionally faceted by user-specified variables and (b) a 2D scatterplot, such as a UMAP or tSNE plot, showing cells colored by the attractors to which they are assigned.

        Parameters
        ----------

        plottingData: pandas DataFrame
            Usually a 2-column dataframe where rows = cells (indices must contain the cells in the training dataset) and columns = UMAP or tSNE (or similar) dimensions
        distanceDF: pandas DataFrame
            Output of assignAttractors, usually. A pandas DataFrame where rows = cells (indices must contain the cells in the training dataset) and at least one column (specified by attractorsAttribute parameter, see below) contains the attractors to which the cells have been assigned
        attractorsAttribute: str
            The column in distanceDF which contains the attractors to which the cells have been assigned. Default is 'decider', as this is the name in the output of assignAttractors
        allAttractors: bool
            True if all assigned attractors are to be included in the frequency and UMAP plots. Should usually be 'True' for making preliminary analysis plots, and 'False' while making publication-quality plots, for ease of visualization and to leave out infrequent attractors.
        numberOfAttractorsToShow: int
            The number of attractors to be shown in the analysis plots. Used when allAttractors is False. Only the top (numberOfAttractorsToShow) most frequent attractors will be shown. The others will be collapsed into a single category "Others".
        cmap: str
            The matplotlib color palette to be used for the analysis plots. Consider using 'colorblind' if the number of categories/attractors is small and a continuous palette such as 'Blues' or 'Greens' if the number is large. See matplotlib documentation for more options on color palettes.
        makeFrequencyPlots: bool
            Should frequency plots be made?
        frequencyGrouping: str
            The variable in the plottingData columns which should be used for faceting the frequency plots. Eg. sample, batch, disease/control, etc.
        freqplotsFile: str
            Path to output PDF file for frequency plots. Default plots to the standard output.
        makeUMAP: bool
            Should UMAP plots be made?
        umapFile: str
            Path to output PDF file for UMAP plots. Default plots to the standard output.

        Returns
        -------
            None
        """

        if makeFrequencyPlots:
            # Make attractor frequency plots
            sns.set_context("talk")
            plt.figure(figsize=(5, 5))
            plottingData.loc[:, "Attractors"] = distanceDF[attractorsAttribute].astype(
                "category"
            )
            plottingData.loc[:, "Attractors"].cat.categories = [
                g for g in range(len(plottingData.loc[:, "Attractors"].cat.categories))
            ]
            if isinstance(allAttractors, list):
                plottingData.loc[:, "Attractors"] = [
                    temp if temp in allAttractors else "Other"
                    for temp in plottingData.loc[:, "Attractors"]
                ]
                palette = {}
                attrColors = {
                    key: value
                    for key, value in zip(
                        allAttractors,
                        list(
                            sns.color_palette(
                                "colorblind", n_colors=len(set(allAttractors))
                            )
                        ),
                    )
                }
                for attr in set(plottingData.Attractors):
                    if attr in allAttractors:
                        palette[attr] = attrColors[attr]
                    else:
                        palette[attr] = (0.75, 0.75, 0.75)
            else:
                plottingData.Attractors = [
                    temp
                    if temp
                    in plottingData.loc[:, "Attractors"]
                    .value_counts()[: numberOfAttractorsToShow + 1]
                    .index.tolist()
                    else "Other"
                    for temp in plottingData.loc[:, "Attractors"]
                ]
                allAttractors = (
                    plottingData.loc[:, "Attractors"]
                    .value_counts()[: numberOfAttractorsToShow + 1]
                    .index.tolist()
                )
                print(allAttractors)
                palette = {}
                attrColors = {
                    key: value
                    for key, value in zip(
                        allAttractors,
                        list(sns.color_palette(cmap, n_colors=len(set(allAttractors)))),
                    )
                }
                for attr in set(plottingData.loc[:, "Attractors"]):
                    if attr in allAttractors:
                        palette[attr] = attrColors[attr]
                    else:
                        palette[attr] = (0.75, 0.75, 0.75)
            if frequencyGrouping != "":
                histData = (
                    plottingData.groupby(by=frequencyGrouping)["Attractors"]
                    .value_counts(normalize=True)
                    .mul(100)
                    .rename("Percent")
                    .reset_index()
                )
                histData.columns = [frequencyGrouping, "Attractors", "Percent"]
                ax = sns.catplot(
                    data=histData,
                    y="Attractors",
                    x="Percent",
                    col=frequencyGrouping,
                    # palette=palette,
                    color="black",
                    dodge=True,
                    orient="h",
                    legend=True,
                    kind="bar",
                    sharey=True,
                    sharex=True,
                )
                ax.set(
                    xlabel="Percentage of Cells\nMapped to Attractors",
                    ylabel="Subject Identifier",
                    xlim=(0, 100),
                )
            else:
                histData = (
                    plottingData.loc[:, "Attractors"]
                    .value_counts(normalize=True)
                    .mul(100)
                    .reset_index()
                )
                histData.columns = ["Attractors", "Percent"]
                ax = sns.catplot(
                    data=histData,
                    y=histData.columns[0],  #'Attractors',
                    x=histData.columns[1],  #'Percent',
                    # palette=palette,
                    color="black",
                    dodge=True,
                    orient="h",
                    legend=True,
                    kind="bar",
                )
                ax.set(
                    xlabel="Percentage of Cells\nMapped to Attractors",
                    ylabel="Attractor Identifier",
                    xlim=(0, 100),
                )
            if freqplotsFile == "":
                plt.show()
                plt.clf()
            else:
                ax.figure.savefig(freqplotsFile, bbox_inches="tight", pad_inches=0.01)
                plt.show()
                plt.clf()
        else:
            print("No attractor frequency plots generated")
        if makeUMAP:
            ##Make attractor UMAP plots
            plt.figure(figsize=(5, 5))
            if frequencyGrouping != "":
                col = frequencyGrouping
            else:
                col = None
            ax = sns.relplot(
                data=plottingData,
                edgecolor="black",
                linewidth=1,
                x=plottingData.columns[0],
                y=plottingData.columns[1],
                hue="Attractors",
                col=col,
                aspect=1,
                palette=palette,
                legend="full",
                s=50,
            )
            # plt.legend(
            #    loc="best", borderaxespad=0.0, ncol=10 #, bbox_to_anchor=(0.5, -0.15),
            # )
            plt.gca().set_aspect("equal")
            # plt.tight_layout()
            if umapFile != "":
                plt.show()
                plt.clf()
            else:
                # get handles and labels for reuse
                label_params = ax.get_legend_handles_labels()
                ax.get_legend().remove()
                ax.figure.savefig(umapFile, bbox_inches="tight", pad_inches=0.01)
                plt.clf()
                # save legend separately
                figl, axl = plt.subplots()
                axl.axis("off")
                axl.legend(*label_params, loc="best", ncol=10, borderaxespad=0)
                figl.savefig(
                    umapFile + "_legend.pdf", bbox_inches="tight", pad_inches=0.01
                )
                plt.clf()
        else:
            print("No UMAP plots generated")

        return None

    def makeBubblePlots(
        pvalues=pd.DataFrame(),
        adjPValueThreshold=0.05,
        wrap=25,
        height=8,
        width=10,
        palette="colorblind",
        saveAsPDF=True,
        outputFile="",
    ):

        """Uses a dataframe of pathway names and adjusted p-values (usually the output of scBONITA pathway analysis) to generate a bubbleplot showing the pathways and the -log10 of their adjusted p-values as calculated by scBONITA

        Parameters
        ----------

        pvalues: pandas DataFrame
            Usually the output of scBONITA pathway analysis A dataframe containing at least the following columns: "Adjusted P value", "Pathway Name", "log10pvals".
        adjPValueThreshold: float
            the adjusted p-value threshold below which dysregulated pathways are shown on the bubbleplot
        wrap: int
            wrap pathway names at this value on the y axis of the bubbleplot
        height: float
            height of image in inches
        width: float
            width of image in inches
        palette: str
            The matplotlib color palette to be used for the analysis plots. Consider using 'colorblind'. See matplotlib
            documentation for more options on color palettes.
        saveAsPDF: bool
            whether to save the bubbleplot as a PDF
        outputFile:
            name of output image file (optional)

        Returns
        -------
            None
        """
        pvalues.loc[:, "log10pvals"] = [
            -1 * np.log10(temp) for temp in pvalues.loc[:, "Adjusted P value"]
        ]
        pvalues = pvalues.loc[pvalues.loc[:, "Adjusted P value"].lt(adjPValueThreshold)]
        pvalues = pvalues.sort_values("Adjusted P value", ascending=0)
        pvalues["Pathway Name"] = pvalues["Pathway Name"].str.wrap(wrap)
        ax = sns.scatterplot(
            data=pvalues,
            s=200,
            x="log10pvals",
            y="Pathway Name",
            alpha=1,
            palette=palette,
        )
        plt.ylabel("Pathway Names")
        plt.xlabel("-log10 (adjusted p-value)")
        axes = plt.gca()
        axes.yaxis.grid(color="grey", linestyle=(0, (5, 10)), linewidth=0.5)
        if saveAsPDF:
            plt.savefig(outputFile, bbox_inches="tight", height=height, width=width)
            plt.show()
            plt.clf()
        else:
            plt.show()
            plt.clf()
        return None

    def ProcessLocalSearchOutput(self, graph, graphName):
        """
        Processes the output of the local search component of scBONITA and returns a text file containing plain-text Boolean rules in the ERS, and other information such as in-degree and out-degree for each node in the network.

        Parameters:
            graphName: name of the graphml file originally provided for rule inference

        Returns:
            A Pandas DataFrame with the following columns:
            'Node': name of the node in the network
            'In_degree': in-degree in the processed network
            'Out_degree': out-degree in the processed network
            'Equivs_len'
            'Graph': the provided graphName
            'Number_of_Nodes': number of nodes in 'graph'
            'importanceScore': the importance score of each node in 'graph' as calculated by scBONITA
            'AvgLocalError': the average error of the local search in scBONITA-RD.
            'plainEquivs': the rules/equivalent rule set from the local search, in human-readable format
            'bitstringEquivs': the rules/equivalent rule set from the local search, in bitstring format - this is the internal representation of rules in scBONITA and may be ignored by users
            'minimalRules': the minimal rule set (minimizing the number of )

        """

        equivs = graphName + "_processed.graphml_equivs1.pickle"
        print(Path(graphName), Path(equivs))
        if (
            Path(graphName).is_file()
            and Path(equivs).is_file()
            and Path(graphName + "_processed.graphml_importanceScores.csv").is_file()
        ):

            errorsName = glob.glob(str(graphName) + "*_localErrors1.pickle")[0]
            localErrors = pickle.load(open(errorsName, "rb"))

            self._singleCell__inherit(
                graph, removeSelfEdges=False, restrictIncomingEdges=False
            )
            equivs = pickle.load(
                open(graphName + "_processed.graphml_equivs1.pickle", "rb")
            )

            equivs = [tuple(equivs[i]) for i in range(0, len(equivs))]
            # print(equivs)
            equivs_len = [len(equivs[i]) for i in range(0, len(equivs))]

            nodeList = list(self.nodeList)
            minimalRules = self.processERS_minimalRule(
                graphName + "_processed.graphml_equivs1.pickle"
            )
            plainEquivs = {}
            bitstringEquivs = {}
            minimalRulesNodes = {}
            for node in range(0, len(nodeList)):
                plainRules = []
                ers = equivs[node]
                for rule in ers:
                    plainRules.append(
                        self._ruleMaker__writeNode(
                            self.nodeList.index(self.nodeList[node]), rule, self
                        )
                    )
                plainEquivs[nodeList[node]] = set(plainRules)
                bitstringEquivs[nodeList[node]] = set(
                    tuple(ers[i]) for i in range(0, len(ers))
                )
                minimalRulesNodes[nodeList[node]] = self._ruleMaker__writeNode(
                    node,
                    minimalRules[
                        self.individualParse[node] : self.individualParse[node + 1]
                    ],
                    self,
                )
            in_degree = [graph.in_degree(node) for node in nodeList]
            out_degree = [graph.out_degree(node) for node in nodeList]
            # inOutRatio = [
            #    float(inDeg + 1) / (outDeg + 1)
            #    for inDeg, outDeg in zip(in_degree, out_degree)
            # ]
            # Get importance score
            importanceScoreDF = pd.read_csv(
                graphName + "_processed.graphml_importanceScores.csv"
            )
            importanceScoreDict = {}
            for n in list(importanceScoreDF["Node"]):
                temp = importanceScoreDF[importanceScoreDF.Node.isin({n})]
                importanceScoreDict[n] = list(temp["Importance Score"])[0]
            tempdf = pd.DataFrame(
                list(zip(nodeList, in_degree, out_degree, equivs_len)),
                columns=["Node", "In_degree", "Out_degree", "Equivs_len"],
            )
            tempdf["Graph"] = graphName
            tempdf["Number_of_Nodes"] = str(len(nodeList))
            tempdf["importanceScore"] = tempdf["Node"].map(importanceScoreDict)
            tempdf["AvgLocalError"] = localErrors
            tempdf["plainEquivs"] = tempdf["Node"].map(plainEquivs)
            tempdf["bitstringEquivs"] = tempdf["Node"].map(bitstringEquivs)
            tempdf["minimalRules"] = tempdf["Node"].map(minimalRulesNodes)
            return tempdf
