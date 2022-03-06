# Class definition for singleCell objects
# Lots of additional functions but the important bits are the ruleMaker and singleCell class definitions
# Author: M.G. Palshikar
# Date: 22May2019
# Updated 15 July with pathway analysis function

import os
from sklearn import preprocessing
import scipy.sparse as sparse
import scipy
import scipy.stats
import numpy as np

# from bioservices.kegg import KEGG
import networkx as nx
import matplotlib.pyplot as plt
import random
import pickle
import time
from scipy.stats.stats import pearsonr, spearmanr
import glob
import sklearn.metrics
import seaborn
from statistics import mean, stdev
import re
import requests
import ctypes as ctypes
import itertools as itertool
from deap import base, creator, gp, tools
from deap import algorithms as algo
import math as math
import copy as copy
from random import random, seed, shuffle, randint, sample, choice
from collections import defaultdict, deque
import gc as gc
from itertools import chain
from operator import attrgetter, itemgetter
import subprocess
import argparse as argparse
import sys

# sys.path.insert(0,'/gpfs/fs1/home/mpalshik/scTest')
import pandas as pd
import csv
import collections

# from testSimulation import *

# start_time = time()
def flatten(l):

    # https://stackoverflow.com/a/2158532
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, (str, bytes)):
            yield flatten(el)
        else:
            yield el


def makeSIF(pathway, keggObject):
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
            type1 = res["entries"][[x["id"] for x in res["entries"]].index(Id1)]["type"]
            type2 = res["entries"][[x["id"] for x in res["entries"]].index(Id2)]["type"]
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
            type1 = res["entries"][[x["id"] for x in res["entries"]].index(Id1)]["type"]
            type2 = res["entries"][[x["id"] for x in res["entries"]].index(Id2)]["type"]
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


def sif_to_digraph(pathwayGenes, sif):
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


def sif_to_digraph2(sif):
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
        elif node1 == node2:
            # don't add self-loops
            continue
        else:
            G.add_edge(str(edge[0]), str(edge[2]), signal=int(edge[1]))
    # graph post-processing
    # remove singletons/isolates
    G.remove_nodes_from(list(nx.isolates(G)))
    # To do: remove complexes, remove dependences of a node on complexes that include that node (which is a form of self-loop)
    return G


def getPathwayName(hsaURL):
    fileReg = re.compile("NAME\s+(\w+.*)")
    pathwayFile = requests.get("http://rest.kegg.jp/get/" + hsaURL, stream=True)
    for line in pathwayFile.iter_lines():
        line = line.decode("utf-8")
        result = fileReg.match(line)
        if result:
            return result.group(1)


# writes rules as a network
# writes rules as a network
def Get_expanded_network(rules, equal_sign="*="):
    """
          The code is written by Gang Yang, Department of Physics, Penn State University if not specified.
    Return the expanded network for a given Boolean network model.
    The Boolean network model is a DiGraph object in the output format of form_network().
    The Boolean network model can be generated through form_network function by reading a text file in the Booleannet format.
    The Boolean rules will first be converted to a disjuctive normal form before generating the expanded network.
    Parameters
    ----------
    Gread     : the given Boolean network model
    prefix='n': prefix to encode the node name to avoid one node's name is a part of another node's name
    suffix='n': suffix to encode the node name to avoid one node's name is a part of another node's name
                e.g. node name '1' will become 'n1n' in the returned result
    equal_sign: the equal sign of the rule in the returned result, whose default value follows the Booleannet format
    Returns
    -------
    The expanded network for the given Boolean network model.
    """
    composite_nodes = []
    G_expand = nx.DiGraph()
    for line in rules:
        child, update_rule = line.split(equal_sign)  # correctly annootate child, rule
        print(child, update_rule)
        if child[0] == "~":  # figure out parity of node
            normal_child = child[1:].strip()
        else:
            normal_child = child[:].strip()
        if "or" in update_rule:
            parents = update_rule.split(" or ")
        else:
            parents = [update_rule]
        parents.sort()
        for parent in parents:
            parent = (
                parent.replace("not ", "~").replace("(", "").replace(")", "").strip()
            )
            if " and " in parent:
                composite_node = parent.replace(" and ", "_").strip()
                composite_nodes.append(composite_node)
                G_expand.add_edge(composite_node, child)
                for component in composite_node.split("_"):
                    G_expand.add_edge(component.strip(), composite_node)
            elif not parent == child:
                G_expand.add_edge(parent, child)
    for node in list(G_expand.nodes()):
        if node[0] == "~" and not "_" in node:
            G_expand.add_edge(node[1:], node)
    nx.set_node_attributes(
        G_expand,
        name="Display Name",
        values={k: " " if k in composite_nodes else k for k in list(G_expand.nodes())},
    )
    nx.set_node_attributes(
        G_expand,
        name="andNode",
        values={k: 1 if k in composite_nodes else 0 for k in list(G_expand.nodes())},
    )

    edgedict = {}
    for edge in list(G_expand.edges()):
        edgedict[edge] = "a"
    nx.set_edge_attributes(G_expand, name="signal", values=edgedict)

    for node in list(G_expand.nodes()):
        if node[0] == "~" and not "_" in node:
            for downstream in G_expand.successors(node):
                G_expand.add_edge(node[1:], downstream, signal="i")
            G_expand.remove_node(node)
    return G_expand.copy()


class ruleMaker:
    def makeToolBox(self, graph):
        # sets up GA toolbox from deap
        # def buildToolbox( individualLength, bitFlipProb, model, params):
        # self.nodeList=list(graph.nodes())
        # self.nodeList.sort()
        weightTup = (-1.0,)  # specify weights of the errors
        for i in range(len(self.nodeList) - 1):
            weightTup += (-1.0,)

        # MAKE TYPES
        creator.create(
            "FitnessMin", base.Fitness, weights=weightTup
        )  # make a fitness minimization function #the objective function has to be MINIMIZED
        creator.create(
            "individual", list, fitness=creator.FitnessMin
        )  # create a class of individuals that are lists of floats

        # INITIALIZATION
        # register our bitsring generator and how to create an individual, population
        toolbox = base.Toolbox()  # build baseline toolbox
        toolbox.register("genRandomBitString", self.genBits)  # , model=self)
        toolbox.register(
            "individual",
            tools.initIterate,
            creator.individual,
            toolbox.genRandomBitString,
        )
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)

        # REGISTER STATISTICS
        # create statistics toolbox and give it functions
        stats = tools.Statistics(key=lambda ind: ind.fitness.values)
        stats.register("avg", np.mean)
        stats.register("std", np.std)
        stats.register("min", np.min)
        stats.register("max", np.max)

        # REGISTER CROSSOVER, MUTATION, AND SELECTION FUNCTIONS
        # finish registering the toolbox functions
        toolbox.register("mate", tools.cxTwoPoint)
        toolbox.register("mutate", tools.mutFlipBit, indpb=self.params.bitFlipProb)
        toolbox.register("select", self.selNSGA2)
        toolbox.register("similar", np.array_equal)

        # ADD TOOLBOX TO OBJECT
        self.toolbox = toolbox
        self.stats = stats

    def getPreds(self, geneList, graph):
        nodes = list(graph)
        # print(nodes)
        corHist = []
        for n in nodes:
            corDict = {}
            n_index = geneList.index(n)
            preds = list(graph.predecessors(n))
            for p in preds:
                p_index = geneList.index(p)
                # corDict[p]=pearsonr(scObj.binMat[n_index,:].todense().tolist()[0], scObj.binMat[p_index,:].todense().tolist()[0])
                # corDict[p]=spearmanr(scObj.binMat[n_index,:].todense().tolist()[0], scObj.binMat[p_index,:].todense().tolist()[0])
                # if(corDict[p][1]<=0.05):
                #     print(n, p, corDict[p])
                # corDict[p]=sklearn.metrics.adjusted_mutual_info_score(scObj.binMat[n_index,:].todense().tolist()[0], scObj.binMat[p_index,:].todense().tolist()[0], average_method="warn")
                corDict[p] = sklearn.metrics.normalized_mutual_info_score(
                    self.binMat[n_index, :].todense().tolist()[0],
                    self.binMat[p_index, :].todense().tolist()[0],
                    average_method="arithmetic",
                )
                # print(n, p, corDict[p])
                corHist.append(corDict[p])
        self.corHist = corHist

    def __init__(
        self,
        graph,
        removeSelfEdges=False,
        restrictIncomingEdges=True,
        maxIncomingEdges=3,
        groundTruth=False,
        graphName="",
    ):
        if maxIncomingEdges < 3:
            print(
                "The maximum number of incoming edges has been set to less than 3. Meaningless results await you."
            )

        nodeList = list(
            graph.nodes
        )  # define the node list simply as the nodes in the graph.

        ruleGraph = nx.empty_graph(0, create_using=nx.DiGraph)  # Create an empty graph
        # remove self loops from the graph
        if removeSelfEdges:
            for node in nodeList:
                repeat = True
                while repeat:
                    repeat = False
                    if node in list(graph.successors(node)):
                        graph.remove_edge(node, node)
                        repeat = True

        self.nodePositions = [
            self.geneList.index(node) for node in nodeList
        ]  # node positions in geneList
        self.nodeList = nodeList
        print("Nodelist: " + str(self.nodeList))
        self.permList = []
        # set up empty lists and dicts for later
        self.rvalues = []  # stores the correlations
        individualParse = (
            []
        )  # list of the number of shadow and nodes that contribute to each node, in order by index num
        andNodeList = (
            []
        )  # a list of the shadow nodes that represent and relations between incoming edge
        andNodeInvertList = (
            []
        )  # keeps track of which incoming nodes for each node need to be inverted
        andLenList = (
            []
        )  # keeps track of how many nodes are coming into each shadow AND node
        nodeDict = (
            {}
        )  # i dentifies names of nodes with their index in the node list- provide name, get index
        possibilityLister = []
        possibilityInverter = []
        succnum = []

        for i in range(0, len(nodeList)):
            nodeDict[
                nodeList[i]
            ] = i  # constructs the node dict so we can easily look up nodes

        counter = int(0)  # keeps track of where we are in the generic individual
        for i in range(0, len(nodeList)):
            predecessors_temp = list(
                graph.predecessors(nodeList[i])
            )  # get NAMES of predecessors of node as documented in the original graph
            successors_temp = list(
                graph.successors(nodeList[i])
            )  # get NAMES of successors of node as documented in the original graph
            succnum.append(len(successors_temp))
            possibilitytemp = [nodeDict[predder] for predder in predecessors_temp]
            possibilityLister.append(list(possibilitytemp))
            # Find correlation between the predecessors and the node
            nodeData = (
                self.binMat[self.nodePositions[i], :].todense().tolist()[0]
            )  # find binarized expression data for node "i"
            predCorr_temp = (
                []
            )  # temporarily store correlations between node "i" and all its predecessors
            for k in predecessors_temp:
                predIndex = self.geneList.index(
                    k
                )  # find index of predecessor in the geneList from the data
                predData = (
                    self.binMat[predIndex, :].todense().tolist()[0]
                )  # find binarized expression data for predecessor "k"
                mi, pvalue = spearmanr(nodeData, predData)
                if np.isnan(mi):
                    predCorr_temp.append(0)
                else:
                    predCorr_temp.append(mi)  # store the calculated correlation
            predecessors_final = sorted(
                zip(predecessors_temp, predCorr_temp),
                reverse=True,
                key=lambda corrs: corrs[1],
            )[
                :3
            ]  # find the top predecessors of the node "i"
            self.rvalues.append(
                sorted(predCorr_temp, reverse=True)[:3]
            )  # stores the correlations
            self.permList.append([pred[0] for pred in predecessors_final])
            for parent in predecessors_final:
                if "interaction" in list(graph[parent[0]][nodeList[i]].keys()):
                    ruleGraph.add_edge(
                        parent[0],
                        nodeList[i],
                        weight=parent[1],
                        activity=graph[parent[0]][nodeList[i]]["interaction"],
                    )
                if "signal" in list(graph[parent[0]][nodeList[i]].keys()):
                    ruleGraph.add_edge(
                        parent[0],
                        nodeList[i],
                        weight=parent[1],
                        activity=graph[parent[0]][nodeList[i]]["signal"],
                    )
            # the following section constructs a list of possible node orders
            # this is accomplished by finding all possible subsets of the list of predecessor nodes

            withNones = zip(
                [nodeList.index(corr_tuple[0]) for corr_tuple in predecessors_final],
                itertool.repeat("empty"),
            )
            possibilities = list(itertool.product(*withNones))
            for j in range(0, len(possibilities)):
                possibilities[j] = list(possibilities[j])
                while "empty" in possibilities[j]:
                    possibilities[j].remove("empty")
                while [] in possibilities[j]:
                    possibilities[j].remove([])
            while [] in possibilities:
                possibilities.remove([])

            # create a list of the activities of each node and store alongside the contributors to each and node for easy reference later
            activities = []  # list to store activities of nodes (a vs i)
            activity = []
            for sequence in possibilities:
                activity = []
                for node in sequence:
                    # check the 'interaction' edge attribute
                    if "interaction" in list(graph[nodeList[node]][nodeList[i]].keys()):
                        if graph[nodeList[node]][nodeList[i]]["interaction"] == "a":
                            activity.append(False)
                        else:
                            if graph[nodeList[node]][nodeList[i]]["interaction"] == "i":
                                activity.append(True)
                            else:
                                if (
                                    graph[nodeList[node]][nodeList[i]]["interaction"]
                                    == "u"
                                ):
                                    print(
                                        "Unknown interaction type, assigning activation..."
                                    )
                                    activity.append(False)
                                else:
                                    if (
                                        graph[nodeList[node]][nodeList[i]][
                                            "interaction"
                                        ]
                                        == "g"
                                    ):
                                        print(
                                            "Group edge/interaction type, assigning activation..."
                                        )
                                        activity.append(False)
                                    else:
                                        print(
                                            "Unknown interaction, assigning activation..."
                                        )
                                        activity.append(False)
                    # check the 'signal' edge attribute
                    if "signal" in list(graph[nodeList[node]][nodeList[i]].keys()):
                        if graph[nodeList[node]][nodeList[i]]["signal"] == "a":
                            activity.append(False)
                        else:
                            if graph[nodeList[node]][nodeList[i]]["signal"] == "i":
                                activity.append(True)
                            else:
                                if graph[nodeList[node]][nodeList[i]]["signal"] == "u":
                                    print(
                                        "Unknown interaction type, assigning activation..."
                                    )
                                    activity.append(False)
                                else:
                                    if (
                                        graph[nodeList[node]][nodeList[i]]["signal"]
                                        == "g"
                                    ):
                                        print(
                                            "Group edge/interaction type, assigning activation..."
                                        )
                                        activity.append(False)
                                    else:
                                        print(
                                            "Unknown interaction, assigning activation..."
                                        )
                                        activity.append(False)
                    # If neither edge attribute is present, assign activation
                    if not "interaction" in list(
                        graph[nodeList[node]][nodeList[i]].keys()
                    ) and not "signal" in list(
                        graph[nodeList[node]][nodeList[i]].keys()
                    ):
                        print("Group edge/interaction type, assigning activation...")
                        activity.append(False)

                activities.append(activity)
            andNodeList.append(possibilities)
            andNodeInvertList.append(activities)
            andLenList.append(len(possibilities))
            possibilityInverter.append(list(activity))
            # construct the list of lengths of possibilties for each node, add to the counter that keeps track of how many bits are necessary
            individualParse.append(counter)
            counter = counter + len(possibilities)

            # Test print statements
            # print("Node: " + self.nodeList[i] + "\nPreds_graph: " + str(sorted(zip(predecessors_temp, predCorr_temp), reverse=True, key=lambda corrs: corrs[1])) + "\nPreds_calculated:  " + str(predecessors_final) + "\nPossibilities: " + str(possibilities)+ "\nActivities/AndNodeInvertList: " + str(activities) + "\nAndNodeList/Possibilities: " + str(possibilities) + "\n****")

        self.size = counter
        individualParse.append(counter)
        self.individualParse = (
            individualParse  # index of start value of current node on the individual
        )
        self.andNodeList = andNodeList  # shadow and node inputs
        self.andNodeInvertList = andNodeInvertList  # keeps track of which incoming nodes for each node need to be inverted
        self.andLenList = (
            andLenList  # keeps track of length of above inputOrderList for each node
        )
        self.possibilityList = possibilityLister  ############
        self.possibilityInverter = possibilityInverter  ##########
        #############
        self.nodeNum = len(nodeList)
        self.params = self.Params()
        self.params.simParams()
        self.makeToolBox(graph)
        self.ruleGraph = ruleGraph
        self.nodeDict = nodeDict  # identifies names of nodes with their index in the node list.. provide name, get index
        self.successorNums = succnum
        # nx.write_graphml(ruleGraph, graphName+"_ruleGraph.graphml")
        print("\nIndividual parse: " + str(self.individualParse))
        print("\nNodelist: " + str(self.nodeList))
        print("\nNode positions: " + str(self.nodePositions))
        print("\nPossibilityList: " + str(self.possibilityList))

    def update_upstream(self, node, newUpstreams):
        withNones = zip(newUpstreams, itertool.repeat("empty"))
        possibilities = list(itertool.product(*withNones))
        for j in range(0, len(possibilities)):
            possibilities[j] = list(possibilities[j])
            while "empty" in possibilities[j]:
                possibilities[j].remove("empty")
            while [] in possibilities[j]:
                possibilities[j].remove([])
        while [] in possibilities:
            possibilities.remove([])
        # create a list of the activities of each node and store alongside the contributors to each and node for easy reference later
        activities = []  # list to store activities of nodes (a vs i)
        for sequence in possibilities:
            activity = []
            for node1 in sequence:
                if (
                    self.possibilityInverter[self.possibilityList[node].index(node1)]
                    == "a"
                ):
                    activity.append(False)
                else:
                    activity.append(True)
            activities.append(activity)
        self.andNodeList[node] = possibilities
        self.andNodeInvertList[node] = activities

    # set up C pointers with correct lengths to pass to simulation software in C
    def updateCpointers(self):
        tempandnoder = []
        tempandinverter = []
        for currentNode in range(len(self.nodeList)):
            tempAndNodes = []
            tempandNodeInvertList = []
            if currentNode < len(self.nodeList):
                tempAndNodes = [
                    xi + [-1] * (3 - len(xi)) for xi in self.andNodeList[currentNode]
                ]
                tempandNodeInvertList = [
                    xi + [-1] * (3 - len(xi))
                    for xi in self.andNodeInvertList[currentNode]
                ]
            while len(tempAndNodes) < 7:
                tempAndNodes.append([0, 0, 0])
                tempandNodeInvertList.append([0, 0, 0])
            tempandnoder.append(tempAndNodes)
            tempandinverter.append(tempandNodeInvertList)
        self.andNodeInvert = np.array(tempandinverter, dtype=np.intc, order="C")
        self.andNodes = np.array(tempandnoder, dtype=np.intc, order="C")

    def genRandBits(self):  # makes a random bitstring
        arr = np.random.randint(2, size=(int(self.size),))
        return list(arr)

    def findEnd(self, node):
        if node == len(self.nodeList) - 1:
            end = self.size
        else:
            end = self.individualParse[node + 1]
        return end

    # executes two point crossover at node junctions
    def cxTwoPointNode(self, ind1, ind2):
        # copied and modified from deap
        # needed to account for bistring only being one of two components of individual
        """Executes a two-point crossover on the input :term:`sequence`
        individuals. The two individuals are modified in place and both keep
        their original length.
        :returns: A tuple of two individuals.
        This function uses the :func:`~random.randint` function from the Python
        base :mod:`random` module.

        Modified to cross over between rules
        """
        size = len(ind1[0].nodeList)
        cxpointer1 = randint(1, size)
        cxpointer2 = randint(1, size - 1)
        # make sure pointers are in right order
        if cxpointer2 >= cxpointer1:
            cxpointer2 += 1
        else:  # Swap the two cx points
            cxpointer1, cxpointer2 = cxpointer2, cxpointer1
        cxpoint1 = ind1[0].individualParse[cxpointer1]
        cxpoint2 = ind1[0].individualParse[cxpointer2]
        # cross over both bitlists and the andNodeLists (as well as andNodeInvertLists)
        ind1[1][cxpoint1:cxpoint2], ind2[1][cxpoint1:cxpoint2] = (
            ind2[1][cxpoint1:cxpoint2],
            ind1[1][cxpoint1:cxpoint2],
        )
        (
            ind1[0].andNodeList[cxpointer1:cxpointer2],
            ind2[0].andNodeList[cxpointer1:cxpointer2],
        ) = (
            ind2[0].andNodeList[cxpointer1:cxpointer2],
            ind1[0].andNodeList[cxpointer1:cxpointer2],
        )
        (
            ind1[0].andNodeInvertList[cxpointer1:cxpointer2],
            ind2[0].andNodeInvertList[cxpointer1:cxpointer2],
        ) = (
            ind2[0].andNodeInvertList[cxpointer1:cxpointer2],
            ind1[0].andNodeInvertList[cxpointer1:cxpointer2],
        )
        # update the arrays seen by C code updateBool
        ind1[0].updateCpointers()
        ind2[0].updateCpointers()
        return ind1, ind2

    # finds the lowest error individual in a population
    def findPopBest(self, population):
        saveVal = -1
        minny = float("Inf")
        for i in range(len(population)):
            if np.sum(population[i].fitness.values) < minny:
                minny = np.sum(population[i].fitness.values)
                saveVal = i
        ultimate = population[saveVal]
        minvals = population[saveVal].fitness.values
        return minvals, ultimate[1], ultimate[0]

    # init value generator for EBNs
    def genEBNInitValues(self, individual, model, sampleProbs):
        # return [True if (random()<sampleProbs[node]) else False for node in range(0,len(sampleProbs))]
        initValues = np.zeros(len(self.sampleList), dtype=np.intc, order="C")
        # for node in range(0,len(sampleProbs)):
        #    if random()<sampleProbs[node]:
        #        initValues[node]=1
        for node in range(0, len(sampleProbs)):
            initValues[node] = sampleProbs[node]
        return initValues

    # NP simulation code for synchronous simulation... fast bc in C
    def NP(self, individual, model, cells, params, KOs, KIs, scSyncBoolC):

        cellArray = []

        # set up knockin and knockout lists
        # knockins=np.array(np.zeros(len(model.nodeList)),dtype=np.intc, order='C')
        # knockouts=np.array(np.zeros(len(model.nodeList)),dtype=np.intc, order='C')

        knockins = np.zeros(len(model.nodeList), dtype=np.intc, order="C")
        knockouts = np.zeros(len(model.nodeList), dtype=np.intc, order="C")

        for knocker in KOs:
            knockouts[knocker] = 1
        for knocker in KIs:
            knockins[knocker] = 1

        # put objects in correct format for passing to C
        nodeIndividual = np.array(individual, dtype=np.intc, order="C")
        indLen = len(nodeIndividual)
        andNodes = np.array(model.andNodes, dtype=np.intc, order="C")
        nodeNum = len(model.nodeList)
        andNodeInvert = np.array(model.andNodeInvert, dtype=np.intc, order="C")
        # print("nodenum: ", nodeNum)
        individualParse = np.array(model.individualParse, dtype=np.intc, order="C")
        andLenList = np.array(model.andLenList, dtype=np.intc, order="C")
        nodePositions1 = model.nodePositions
        nodePositionsC = np.array(nodePositions1, dtype=np.intc, order="C")
        # print(nodePositionsC)
        simSteps = self.params.simSteps
        lenSamples1 = len(model.sampleList)
        binMatC1 = self.binMat.toarray(order="C")
        # binMatC2=np.transpose(binMatC1)
        # binMatC3=np.array(binMatC2, order="C",dtype=np.intc)
        binMatC3 = np.transpose(
            np.array(copy.deepcopy(binMatC1), order="C", dtype=np.intc)
        )
        # print(binMatC3)
        # print(binMatC3.shape)
        binMatCPointer = ctypes.c_void_p(
            binMatC3.ctypes.data
        )  # put input array as C pointer

        # convert objects into C pointers
        nodeIndividual1 = ctypes.c_void_p(nodeIndividual.ctypes.data)
        indLen1 = ctypes.c_void_p(indLen)
        andNodes1 = ctypes.c_void_p(andNodes.ctypes.data)
        individualParse1 = ctypes.c_void_p(individualParse.ctypes.data)
        andLenList1 = ctypes.c_void_p(andLenList.ctypes.data)
        andNodeInvertList1 = ctypes.c_void_p(andNodeInvert.ctypes.data)
        nodeNum1 = ctypes.c_void_p(nodeNum)

        simSteps1 = ctypes.c_void_p(simSteps)
        knockouts1 = ctypes.c_void_p(knockouts.ctypes.data)
        knockins1 = ctypes.c_void_p(knockins.ctypes.data)

        nodePositionsCPointer = ctypes.c_void_p(nodePositionsC.ctypes.data)

        vals = np.full(
            shape=(self.maxNodes, self.params.simSteps, self.maxSamples),
            fill_value=2,
            dtype=np.intc,
            order="C",
        )  # initiate output array - sim data is nodes * sim steps * cells. Max sim steps hard coded to 200
        valsubmit = ctypes.c_void_p(vals.ctypes.data)  # put output array into C pointer
        lenSamples = ctypes.c_void_p(lenSamples1)
        # print(nodeIndividual)
        # print(andNodes)
        # print(andNodeInvert)
        # print("lensamples: ")
        # print(lenSamples)

        scSyncBoolC(
            valsubmit,
            nodeIndividual1,
            indLen1,
            nodeNum1,
            andLenList1,
            individualParse1,
            andNodes1,
            andNodeInvertList1,
            simSteps1,
            knockouts1,
            knockins1,
            lenSamples,
            binMatCPointer,
            nodePositionsCPointer,
        )
        # print("Simulation trajectory: ")
        # print(vals.shape)
        # print(vals[:, :, 0])

        # vals2=copy.deepcopy(vals)
        # print(vals2)
        # return [( (1.*np.sum(col))/cells) for col in zip(*cellArray)]
        # return [((1.*np.sum(col))/cells) for col in vals2]
        return vals

    # calculate  fitness for an individual

    def evaluateByNode(
        self,
        individual,
        cells,
        model,
        sss,
        params,
        KOlist,
        KIlist,
        scSyncBoolC,
        localSearch=False,
    ):
        # sss=model.nodePositions
        # boolValues=[self.NP(list(individual), model, cells, model.initValueList[i], self.params, KOlist[i], KIlist[i], boolC) for i in range(len(sss))]
        # boolValues=[self.NP(list(individual), model, cells, model.initValueList[i], self.params, KOlist[i], KIlist[i], scSyncBoolC) for i in range(len(sss))]
        # boolValues=self.NP(individual, scTest, scTest.params.cells, scTest.params,scTest.knockinLists, scTest.knockinLists, scSyncBoolC)
        boolValues = self.NP(
            individual,
            model,
            self.params.cells,
            self.params,
            KOlist,
            KIlist,
            scSyncBoolC,
        )
        # self.nodePositions are node positions in geneList
        # if localSearch:
        #    print("Boolvalues: ")
        #    print(boolValues)
        # print("sss: ")
        # print(sss)
        # return tuple([np.sum([(boolValues[j][i1]-self.binMat[i2, j])**2 for j in range(0,len(sss))]) for i2, i1 in zip(self.nodePositions, range(0, len(model.nodeList)))]) #fixed now
        # Get position of model.nodeList nodes in the geneList - these are the column indices for binMat
        # Get position of sss (ie, samples) in the sampleList - these are the row indices for binMat
        # check structure of boolValues
        return self.calculateErrorFromSimData(
            simData=boolValues, localSearch=localSearch, model=model
        )

    def calculateErrorFromSimData(self, simData, localSearch, model):
        allCellErrors = []
        for node in list(model.nodePositions):
            # print(node)
            simAvg_k = []  # postulated attractor average for the starting state k
            # for simStep in range(self.params.simSteps - 3, self.params.simSteps): #iterate over the last three steps of the simulation trajectory
            for simStep in range(
                self.params.simSteps - 1, self.params.simSteps
            ):  # iterate over the last three steps of the simulation trajectory
                for k in range(
                    0, len(model.sampleList)
                ):  # iterate over all cells, ie, over all starting states
                    # print([k, simStep, node])
                    simAvg_k.append(simData[node, simStep, k])
            # print(simAvg_k)
            # simAvg_k = [float(i)/3 for i in simAvg_k]
            # print(simAvg_k)
            # Get start value of cell from binMat
            # cellInitValue = self.binMat[:, m]
            cellInitValue = self.binMat[node, :]
            cellInitValue = cellInitValue.toarray()
            cellInitValue = cellInitValue.tolist()[0]
            # cellInitValue = [cellTemp[0] for cellTemp in cellInitValue]

            # Calculate error for cell k
            error_k = [
                abs(tempCell - tempAttr)
                for tempCell, tempAttr in zip(cellInitValue, simAvg_k)
            ]
            # print('Error_k: '+str(error_k[1:5]) + '\tcellInitValue: ' + str(cellInitValue[1:5])+ '\tSim_avg_k: ' + str(simAvg_k[1:5]))
            allCellErrors.append(error_k)
        # print('Error_k: '+str(error_k[1:5]) + '\tcellInitValue: ' + str(cellInitValue[1:5])+ '\tSim_avg_k: ' + str(simAvg_k[1:5]))
        # print(allCellErrors[1:5]) #allCellErrors is now a list of lists. Each list corresponds to a cell (or a start state, and each value in that list corresponds to the error for a node)
        if localSearch:
            allNodeErrors = np.array(allCellErrors)
            # print("allNodeErrors:\t"+str(allNodeErrors.shape))
            allNodeErrors = [
                sum(allNodeErrors[i, :]) / len(model.sampleList)
                for i in range(0, len(model.nodePositions))
            ]
            # allNodeErrors=mean(allNodeErrors)
            # print(allNodeErrors[0:5])
            return allNodeErrors

        else:
            allCellErrors = [sum(temp) for temp in allCellErrors]
            allCellErrors = [
                float(temp) / len(model.sampleList) for temp in allCellErrors
            ]
            return allCellErrors

    def calculateErrorFromAttractors():
        attractors = getAttractors_python(vals2D, updateBooler, nodeNum=200)

    def calculateErrorFromSimData2(self, simData, localSearch, model):
        if localSearch:
            finalError = [0.0] * len(model.nodePositions)
            for node in self.nodePositions:
                print(node)
                d_ij = self.binMat.tocsr()[
                    node, range(0, len(model.sampleList))
                ].toarray()[
                    0
                ]  # value of node 'node' in all samples
                o_ij = simData[
                    node, self.params.simSteps - 1, range(0, len(model.sampleList))
                ]  # value of node in the
                print(simData[node, :, :])
                d_ij = list(d_ij)
                o_ij = list(o_ij)
                errorTemp = [
                    abs(tempCell - tempAttr) for tempCell, tempAttr in zip(d_ij, o_ij)
                ]
                # errorTemp=abs(d_ij - o_ij)
                print(
                    "error_k "
                    + str(list(errorTemp))
                    + "\t"
                    + str((d_ij))
                    + "\t"
                    + str((o_ij))
                )
                # print(list(errorTemp))
                if errorTemp[node] != 0:
                    errorTemp = [float(60000)] * len(
                        model.nodePositions
                    )  # [float(10000)] * len(model.nodePositions)
                    finalError[node] = sum(errorTemp)
                else:
                    finalError[node] = sum(errorTemp)
            print(sum(finalError))
            return finalError
        else:
            finalError = [0.0] * len(model.nodePositions)
            for k in range(
                0, len(model.sampleList)
            ):  # iterate over all cells, ie, over all starting states
                error = 0.0
                for node in range(0, len(self.nodePositions)):
                    # print(self.nodeList[node]) #get gene index
                    # print(self.geneList[node])
                    d_ij = self.binMat[node, k]  # value of node 'node' in sample k
                    o_ij = simData[
                        node, self.params.simSteps - 1, k
                    ]  # value of node in the final step of simulation starting from k
                    errorTemp = abs(d_ij - o_ij)
                    print(
                        "error_k "
                        + str(errorTemp)
                        + "\t"
                        + str(d_ij)
                        + "\t"
                        + str(o_ij)
                    )
                    error = errorTemp + error
                error = float(error) / len(self.nodePositions)
                finalError[node] = error
            # print(finalError)
            return finalError

    def calculateImportanceScoreFromSimData(self, simData1, simData2):
        # print(self.binMat.shape)
        allCellErrors = []
        for k in range(
            0, len(self.sampleList)
        ):  # iterate over all cells, ie, over all starting states
            simAvg_k1 = []  # postulated attractor average for the starting state k
            simAvg_k2 = []
            for simStep in range(
                self.params.simSteps - 1, self.params.simSteps
            ):  # iterate over the last three steps of the simulation trajectory
                for node in self.nodePositions:
                    # print([k, simStep, node])
                    simAvg_k1.append(simData1[node, simStep, k])
                    simAvg_k2.append(simData2[node, simStep, k])
            # print(simAvg_k)
            simAvg_k1 = [float(i) for i in simAvg_k1]
            simAvg_k2 = [float(i) for i in simAvg_k2]
            error_k = [
                abs(tempCell1 - tempCell2)
                for tempCell1, tempCell2 in zip(simAvg_k1, simAvg_k2)
            ]
            allCellErrors.append(error_k)
        allCellErrors = [sum(temp) for temp in allCellErrors]
        allCellErrors = [float(temp) / len(self.nodeList) for temp in allCellErrors]
        # print(allCellErrors)
        allCellErrors = sum(allCellErrors)
        return allCellErrors

    # generates list of offspring to be compared... decides to do crossover or mutation
    def varOrAdaptive(
        self, population, toolbox, lambda_, cxpb, mutpb, genfrac, mutModel
    ):
        # def varOrAdaptive(population, toolbox, model, lambda_, cxpb, mutpb, genfrac, mutModel):
        # algorithm for generating a list of offspring... copied and pasted from DEAP with modification for adaptive mutation
        assert (cxpb + mutpb) <= 1.0, (
            "The sum of the crossover and mutation "
            "probabilities must be smaller or equal to 1.0."
        )
        offspring = []
        for _ in range(lambda_):
            op_choice = random()
            if op_choice < cxpb:  # Apply crossover
                inds = []
                for samp in sample(population, 2):
                    ind = toolbox.clone(samp)
                    inds.append(ind)
                ind1, ind2 = inds
                ind1, ind2 = self.cxTwoPointNode(ind1, ind2)
                del ind1.fitness.values
                offspring.append(ind1)
            elif op_choice < cxpb + mutpb:  # Apply mutation
                ind = toolbox.clone(choice(population))
                (ind,) = self.mutFlipBitAdapt(ind, genfrac, mutModel)
                del ind.fitness.values
                offspring.append(ind)
            else:  # shouldn't happen... clone existing individual
                offspring.append(choice(population))
        return offspring

    # select node to mutate
    def selectMutNode(self, errors):
        # print("normerrors: " + str(errors))
        normerrors = [
            1.0 * error / np.sum(errors) for error in errors
        ]  # normalize errors to get a probability that the node  is modified
        probs = np.cumsum(normerrors)
        randy = random()  # randomly select a node to mutate
        return next(i for i in range(len(probs)) if probs[i] > randy)

    # mutation algorithm
    def mutFlipBitAdapt(self, indyIn, genfrac, mutModel):
        errors = list(indyIn.fitness.values)  # get errors
        individual = indyIn[1]
        model = indyIn[0]
        # get rid of errors in nodes that can't be changed
        errorNodes = 0
        for j in range(len(errors)):
            if model.andLenList[j] < 2:
                errors[j] = 0
            else:
                errorNodes = errorNodes + 1

        if np.sum(errors) < 0.05 * errorNodes or errorNodes == 0:
            # condition selection on number of incoming edges + downstream edges
            pseudoerrors = [
                len(model.possibilityList[i])
                if model.successorNums[i] == 0
                else len(model.possibilityList[i]) * model.successorNums[i]
                for i in range(len(model.nodeList))
            ]
            # zero out nodes that can't be changed
            for j in range(len(pseudoerrors)):
                if model.andLenList[j] < 2:
                    pseudoerrors[j] = 0
            focusNode = self.selectMutNode(pseudoerrors)
        else:
            # if errors are relatively high, focus on nodes that fit the worst and have highest in-degree
            # calculate probabilities for mutating each node
            for i in range(len(errors)):
                temper = model.successorNums[i]
                if temper == 0:
                    errors[i] = errors[i] * len(model.possibilityList[i])
                else:
                    errors[i] = errors[i] * len(model.possibilityList[i]) * temper
            focusNode = self.selectMutNode(errors)
        # perform mutation
        if model.andLenList[focusNode] > 1:
            # find ends of the node of interest in the individual
            start = model.individualParse[focusNode]
            end = model.findEnd(focusNode)
            # mutate the inputs some of the time
            if len(model.possibilityList[focusNode]) > 3 and random() < mutModel:
                temppermup = []  # temporary upstream nodes
                upstreamAdders = list(model.possibilityList[focusNode])
                rvals = list(model.rvalues[focusNode])
                while len(temppermup) < 3:
                    randy = random()  # randomly select a node to mutate
                    tempsum = sum(rvals)
                    if tempsum == 0:
                        addNoder = randint(
                            0, len(rvals) - 1
                        )  # int(math.floor(random()*len(upstreamAdders)))
                        print(addNoder)
                    else:
                        recalc = np.cumsum([1.0 * rval / tempsum for rval in rvals])
                        print(recalc)
                        addNoder = next(
                            i for i in range(len(recalc)) if recalc[i] > randy
                        )
                        print(addNoder)
                    temppermup.append(upstreamAdders.pop(addNoder))
                    print(rvals)
                    rvals.pop(addNoder)
                model.update_upstream(focusNode, temppermup)
                model.updateCpointers()
            for i in range(start, end):
                # print("i: " + str(i))
                if random() < 2 / (end - start + 1):
                    individual[i] = 1
                else:
                    individual[i] = 0
            # ensure that there is at least one shadow and node turned on
            if np.sum(individual[start:end]) == 0:
                individual[start] = 1
            indyIn[0] = model
            indyIn[1] = individual
        else:
            print("did not actually check")
        return (indyIn,)

    def genBits(self):
        # generate random bitlist
        startInd = list(self.genRandBits())
        counter = 0
        # make sure bitlist isn't zero
        while np.sum(startInd) == 0 and counter < float("Inf"):
            startInd = list(self.genRandBits())
            counter += 1
        # go through nodes and make sure that there are 1-5 ones in the random list
        for node in range(0, len(self.nodeList)):
            end = self.findEnd(node)
            start = self.individualParse[node]
            if (end - start) > 1:
                counter = 0
                while np.sum(startInd[start:end]) > 5 and counter < float("Inf"):
                    chosen = math.floor(random() * (end - start))
                    startInd[start + int(chosen)] = 0
                    counter += 1
                if np.sum(startInd[start:end]) == 0:
                    chosen = math.floor(random() * (end - start))
                    startInd[start + int(chosen)] = 1
            elif (end - start) == 1:
                startInd[start] = 1
        return [copy.deepcopy(self), startInd]

    # taken from deap and modified slightly to make pareto sorting less strict
    def sortNondominatedAdapt(self, individuals, k, first_front_only=False):
        """Sort the first *k* *individuals* into different nondomination levels
        using the "Fast Nondominated Sorting Approach" proposed by Deb et al.,
        see [Deb2002]_. This algorithm has a time complexity of :math:`O(MN^2)`,
        where :math:`M` is the number of objectives and :math:`N` the number of
        individuals.

        :param individuals: A list of individuals to select from.
        :param k: The number of individuals to select.
        :param first_front_only: If :obj:`True` sort only the first front and
                                exit.
        :returns: A list of Pareto fronts (lists), the first list includes
                nondominated individuals.
        .. [Deb2002] Deb, Pratab, Agarwal, and Meyarivan, "A fast elitist
        non-dominated sorting genetic algorithm for multi-objective
        optimization: NSGA-II", 2002.
        """
        if k == 0:
            return []

        map_fit_ind = defaultdict(list)
        for ind in individuals:
            map_fit_ind[ind.fitness].append(ind)

        fits = list(map_fit_ind)
        current_front = []
        next_front = []
        dominating_fits = defaultdict(int)
        dominated_fits = defaultdict(list)

        # Rank first Pareto front
        for i, fit_i in enumerate(fits):
            for fit_j in fits[i + 1 :]:
                if self.dominated(fit_i, fit_j):
                    dominating_fits[fit_j] += 1
                    dominated_fits[fit_i].append(fit_j)
                elif self.dominated(fit_j, fit_i):
                    dominating_fits[fit_i] += 1
                    dominated_fits[fit_j].append(fit_i)
            if dominating_fits[fit_i] == 0:
                current_front.append(fit_i)

        fronts = [[]]
        for fit in current_front:
            fronts[-1].extend(map_fit_ind[fit])
        pareto_sorted = len(fronts[-1])

        # Rank the next front until all individuals are sorted or
        # the given number of individual are sorted.
        if not first_front_only:
            N = min(len(individuals), k)
            while pareto_sorted < N:
                fronts.append([])
                for fit_p in current_front:
                    for fit_d in dominated_fits[fit_p]:
                        dominating_fits[fit_d] -= 1
                        if dominating_fits[fit_d] == 0:
                            next_front.append(fit_d)
                            pareto_sorted += len(map_fit_ind[fit_d])
                            fronts[-1].extend(map_fit_ind[fit_d])
                current_front = next_front
                next_front = []
        return fronts

    # taken from deap and modified slightly to make pareto sorting less strict
    def dominated(self, ind1, ind2):
        """Return true if each objective of *self* is not strictly worse than
        the corresponding objective of *other* and at least one objective is
        strictly better.
        :param obj: Slice indicating on which objectives the domination is
                    tested. The default value is `slice(None)`, representing
                    every objectives.
        """
        not_equal = False
        mean1 = np.mean(ind1.wvalues)
        mean2 = np.mean(ind2.wvalues)
        std1 = np.std(ind1.wvalues)
        if mean1 > mean2:
            not_equal = True
        elif mean1 < mean2:
            return False
        return not_equal

    # taken from deap
    def assignCrowdingDist(self, individuals):
        """Assign a crowding distance to each individual's fitness. The
        crowding distance can be retrieve via the :attr:`crowding_dist`
        attribute of each individual's fitness.
        """
        if len(individuals) == 0:
            return

        distances = [0.0] * len(individuals)
        crowd = [(ind.fitness.values, i) for i, ind in enumerate(individuals)]

        nobj = len(individuals[0].fitness.values)

        for i in range(nobj):
            crowd.sort(key=lambda element: element[0][i])
            distances[crowd[0][1]] = float("inf")
            distances[crowd[-1][1]] = float("inf")
            if crowd[-1][0][i] == crowd[0][0][i]:
                continue
            norm = nobj * float(crowd[-1][0][i] - crowd[0][0][i])
            for prev, cur, next in zip(crowd[:-2], crowd[1:-1], crowd[2:]):
                distances[cur[1]] += 1.0 * (next[0][i] - prev[0][i]) / norm

        for i, dist in enumerate(distances):
            individuals[i].fitness.crowding_dist = dist

    def selNSGA2(self, individuals, k):
        # NSGA2 selection taken from deap
        """Apply NSGA-II selection operator on the *individuals*. Usually, the
        size of *individuals* will be larger than *k* because any individual
        present in *individuals* will appear in the returned list at most once.
        Having the size of *individuals* equals to *k* will have no effect other
        than sorting the population according to their front rank. The
        list returned contains references to the input *individuals*. For more
        details on the NSGA-II operator see [Deb2002]_.

        :param individuals: A list of individuals to select from.
        :param k: The number of individuals to select.
        :returns: A list of selected individuals.

        .. [Deb2002] Deb, Pratab, Agarwal, and Meyarivan, "A fast elitist
        non-dominated sorting genetic algorithm for multi-objective
        optimization: NSGA-II", 2002.
        """
        pareto_fronts = self.sortNondominatedAdapt(individuals, k)
        for front in pareto_fronts:
            self.assignCrowdingDist(front)

        chosen = list(chain(*pareto_fronts[:-1]))
        k = k - len(chosen)
        if k > 0:
            sorted_front = sorted(
                pareto_fronts[-1], key=attrgetter("fitness.crowding_dist"), reverse=True
            )
            chosen.extend(sorted_front[:k])

        return chosen

    # master GA algorithm
    def eaMuPlusLambdaAdaptive(self, scSyncBoolC, graph, verbose=True):
        params = self.params
        toolbox = self.toolbox
        mutModel = self.params.mutModel
        logbook = tools.Logbook()
        mu = self.params.mu
        lambda_ = self.params.lambd
        stats = self.stats
        cxpb = self.params.crossoverProb
        mutpb = self.params.mutationProb
        ngen = self.params.generations
        sampleList = self.sampleList
        KOlist = self.knockoutLists
        KIlist = self.knockinLists

        # print(mu, lambda_, cxpb, mutpb, ngen, sampleList[0:4], len(KOlist), len(KIlist))

        population = self.toolbox.population(n=self.params.popSize)

        logbook.header = ["gen", "nevals"] + (self.stats.fields if self.stats else [])
        lastcheck = []
        modellist = []
        fitnesslist = []
        popList = []

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in population if not ind.fitness.valid]
        print("Invalid individuals: " + str(invalid_ind))
        updateBooler = ctypes.cdll.LoadLibrary("./simulator.so")
        scSyncBoolC = updateBooler.scSyncBool

        fitnesses = [
            indy[0].evaluateByNode(
                indy[1],
                indy[0].params.cells,
                indy[0],
                indy[0].nodePositions,
                indy[0].params,
                KOlist,
                KIlist,
                scSyncBoolC,
            )
            for indy in invalid_ind
        ]

        print("Fitnesses: " + str(fitnesses))

        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        fitnesslist.append([list(ind.fitness.values) for ind in population])
        popList.append([list(inder[1]) for inder in population])
        modellist.append(
            [
                [
                    (modeler[0].size),
                    list(modeler[0].nodeList),
                    list(modeler[0].individualParse),
                    list(modeler[0].andNodeList),
                    list(modeler[0].andNodeInvertList),
                    list(modeler[0].andLenList),
                    list(modeler[0].nodeList),
                    dict(modeler[0].nodeDict),
                ]
                for modeler in population
            ]
        )

        record = stats.compile(population) if stats is not None else {}
        logbook.record(gen=0, nevals=len(invalid_ind), **record)
        if verbose:
            print(logbook.stream)

        breaker = False
        for ind in population:
            if np.sum(ind.fitness.values) < 0.01 * len(ind.fitness.values):
                breaker = True
        if breaker:
            # outputList=[fitnesslist, popList, modellist]
            # pickle.dump( outputList, open(graph+"_pops.pickle", "wb" ) )
            return population, logbook

        # Begin the generational process
        for gen in range(1, ngen + 1):
            # offspring = self.varOrAdaptive(population, toolbox, lambda_, .5+.5*(1.-1.*gen/ngen), (.5*gen/ngen), (1.*gen/ngen),mutModel)
            # varOrAdaptive(self, population, toolbox, lambda_, cxpb, mutpb, genfrac, mutModel)
            offspring = self.varOrAdaptive(
                population, toolbox, lambda_, cxpb, mutpb, (1.0 * gen / ngen), mutModel
            )
            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            fitnesses = [
                self.evaluateByNode(
                    indy[1],
                    self.params.cells,
                    indy[0],
                    self.nodePositions,
                    self.params,
                    KOlist,
                    KIlist,
                    scSyncBoolC,
                )
                for indy in invalid_ind
            ]
            # fitnesses=[indy[0].evaluateByNode(indy[1], indy[0].params.cells, indy[0], indy[0].nodePositions, indy[0].params, KOlist, KIlist,scSyncBoolC) for indy in invalid_ind]
            print(fitnesses[0][0:5])
            print(invalid_ind[0:5][1])
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit
            # Select the next generation population
            population[:] = toolbox.select(population + offspring, mu)
            fitnesslist.append([list(ind.fitness.values) for ind in population])
            popList.append([list(inder[1]) for inder in population])
            modellist.append(
                [
                    [
                        (modeler[0].size),
                        list(modeler[0].nodeList),
                        list(modeler[0].individualParse),
                        list(modeler[0].andNodeList),
                        list(modeler[0].andNodeInvertList),
                        list(modeler[0].andLenList),
                        list(modeler[0].nodeList),
                        dict(modeler[0].nodeDict),
                    ]
                    for modeler in population
                ]
            )

            # Update the statistics with the new population
            record = stats.compile(population) if stats is not None else {}
            logbook.record(gen=gen, nevals=len(invalid_ind), **record)
            if verbose:
                print(logbook.stream)
            breaker = False
            for ind in population:
                if np.sum(ind.fitness.values) < 0.01 * len(ind.fitness.values):
                    breaker = True
                    saveInd = ind
            if breaker:
                errorTemp = saveInd.fitness.values
                for value in errorTemp:
                    if value > 0.1:
                        breaker = False
            if breaker:
                outputList = [fitnesslist, popList, modellist]
                # pickle.dump( outputList, open(graph+"_pops.pickle", "wb" ) )
                return population, logbook

        # outputList=[fitnesslist, modellist] #popList,
        # pickle.dump( outputList, open(graph+"_pops.pickle", "wb" ) )

        return population, logbook

    def bitList(self, n, x):
        templist = [1 if digit == "1" else 0 for digit in bin(n)[::-1]]
        while len(templist) < x:
            templist.append(0)
        while (len(templist)) > x:
            templist.pop()
        return templist

    # local search function
    def checkNodePossibilities(
        self, node, indy, newSSS, cellNum, model, params, KOlist, KIlist, scSyncBoolC
    ):
        # print(node)
        tol = 0.0  # .01*len(newSSS) # set tolerance for equivalence
        end = model.findEnd(node)  # find end of model for this node
        start = model.individualParse[node]  # find start of model for this node
        truth = list(indy[start:end])
        equivs = [truth]  # add if
        # print(end)
        # print(start)
        if (end - start) == 0:
            # print("Checkpoint")
            return truth, equivs, equivs, 0.0
        indOptions = []
        indErrors = []
        # iterate over possibilities for this node
        for i in range(1, 2 ** (end - start)):
            tempultimate = list(indy)
            tempInd = model.bitList(i, len(truth))  # self.bitList(i, len(truth))
            # print(tempInd)
            tempultimate[start:end] = tempInd  # set rule to one being checked
            currentsumtemp = self.evaluateByNode(
                tempultimate,
                self.params.cells,
                model,
                model.nodePositions,
                self.params,
                KOlist,
                KIlist,
                scSyncBoolC,
                localSearch=True,
            )  # self.evaluateByNode(tempultimate, self.params.cells, self, self.nodePositions, self.params, KOlist, KIlist, scSyncBoolC, localSearch=True) #gives one error per cell
            # print(currentsumtemp)
            currentsum = currentsumtemp[
                node
            ]  # save the error found #subset complete error to
            indOptions.append(tempInd)
            indErrors.append(currentsum)
            gc.collect()
        minny = min(indErrors)
        # print(indErrors)
        # print(minny)
        equivs = []
        # find the minimum error individual
        for i in range(len(indOptions)):
            if indErrors[i] <= minny + tol:  # changed from < to <=
                equivs.append(indOptions[i])
        truth = equivs[0]
        print(
            "Node\t"
            + str(node)
            + "\t"
            + "indErrors\t"
            + str(indErrors)
            + "\t"
            + "Minny\t"
            + str(minny)
            + "\t"
            + "Len_equivs\t"
            + str(len(equivs))
            + "Equivs\t"
            + str(equivs)
        )

        return (truth, equivs, minny, indErrors)

    def writeModel(self, individual, model):
        # iterate over nodes to generate a BooleanNet representation for the entire model
        addString = ""
        # for i in range(0,len(self.nodeList)):
        for i in range(0, len(model.nodePositions)):
            addString = addString + model.writeNode(
                i,
                individual[model.individualParse[i] : model.individualParse[i + 1]],
                model,
            )
            addString = addString + "\n"
        return addString[:-1]

    def findInEdges(self, model, node):
        """find the incoming edges to each 'and' connection for a given node"""
        inEdges = []
        for lister in model.andNodeList[node]:
            tempTup = tuple(lister)
            inEdges.append(set(tempTup))
        return inEdges

    def simplifyRule(self, rule, inEdges):
        """find the simplest form of a rule"""
        for i in range(len(rule)):
            if rule[i] == 1:
                for j in range(len(rule)):
                    if rule[j] == 1 and not i == j:
                        if inEdges[i].issubset(inEdges[j]):
                            rule[j] = 0
        return rule

    def writeNode(self, currentNode, nodeIndividual, model):
        # write out evaluation instructions in BooleanNet format.
        # This follows the exact same code as updateNode (for switch=0), but writes a string instead of actually updating the values of the nodes
        andNodes = model.andNodeList[
            currentNode
        ]  # find the list of shadow and nodes we must compute before computing value of current nodes
        andNodeInvertList = model.andNodeInvertList[
            currentNode
        ]  # find list of lists of whether input nodes need to be inverted (corresponds to inputOrder)
        writenode = (
            "" + model.nodeList[currentNode] + " *= "
        )  # set up the initial string to use to write node
        # print("\nCurrent Node: ")
        # print(currentNode)
        # print(model.nodeList[currentNode])
        # print("\nAnd nodes: ")
        # print(andNodes)
        # print([model.nodeList[an[0]] for an in andNodes])
        # print("\nAnd node invert list: ")
        # print(andNodeInvertList)
        # print("\n Node individual: ")
        # print(nodeIndividual)
        inEdges = self.findInEdges(model, currentNode)
        nodeIndividual = self.simplifyRule(nodeIndividual, inEdges)
        # print("\n Simplified node individual: ")
        # print(nodeIndividual)

        if model.andLenList[currentNode] == 0 or sum(nodeIndividual) == 0:
            # print(writenode + ' ' + model.nodeList[currentNode])
            return (
                writenode + " " + model.nodeList[currentNode]
            )  # if no inputs, maintain value
        elif len(andNodes) == 1:
            # if only one input, then can either affect or not affect the node. so either keep the value or update to the single input's value
            value = ""
            # if only one input, then set to that number
            if andNodeInvertList[0][0] == 0:
                value = value + model.nodeList[andNodes[0][0]]
            else:
                value = value + "not " + model.nodeList[andNodes[0][0]]
            print(writenode + value)
            return writenode + value
        else:
            # update nodes with more than one input
            # first deal with case of simple logic without need of linear regression
            orset = []
            # go through list of possible shadow and nodes to see which ones actually contribute
            for andindex in range(len(nodeIndividual)):
                newval = "("
                if nodeIndividual[andindex] == 1:
                    # if a shadow and contributes, compute its value using its upstream nodes
                    if andNodeInvertList[andindex][0]:
                        newval = newval + "not "
                    newval = newval + self.nodeList[andNodes[andindex][0]]
                    for addnode in range(1, len(andNodes[andindex])):
                        newval = newval + " and "
                        if andNodeInvertList[andindex][addnode]:
                            newval = newval + " not "
                        newval = newval + self.nodeList[andNodes[andindex][addnode]]
                    orset.append(newval + ")")
                # combine the shadow and nodes with or operations
            writenode = writenode + orset.pop()
            for val in orset:
                writenode = writenode + " or " + val
            # print(writenode)
            return writenode

    def writeNode_BoolNet(self, currentNode, nodeIndividual, model):
        # write out evaluation instructions in BoolNet format.
        # This follows the exact same code as updateNode (for switch=0), but writes a string instead of actually updating the values of the nodes
        andNodes = model.andNodeList[
            currentNode
        ]  # find the list of shadow and nodes we must compute before computing value of current nodes
        andNodeInvertList = model.andNodeInvertList[
            currentNode
        ]  # find list of lists of whether input nodes need to be inverted (corresponds to inputOrder)
        writenode = (
            "" + model.nodeList[currentNode] + " , "
        )  # set up the initial string to use to write node
        inEdges = self.findInEdges(model, currentNode)
        nodeIndividual = self.simplifyRule(nodeIndividual, inEdges)
        if model.andLenList[currentNode] == 0 or sum(nodeIndividual) == 0:
            return (
                writenode + " " + model.nodeList[currentNode]
            )  # if no inputs, maintain value
        elif len(andNodes) == 1:
            # if only one input, then can either affect or not affect the node. so either keep the value or update to the single input's value
            value = ""
            # if only one input, then set to that number
            if andNodeInvertList[0][0] == 0:
                value = value + model.nodeList[andNodes[0][0]]
            else:
                value = value + "!" + model.nodeList[andNodes[0][0]]
            print(writenode + value)
            return writenode + value
        else:
            # update nodes with more than one input
            # first deal with case of simple logic without need of linear regression
            orset = []
            # go through list of possible shadow and nodes to see which ones actually contribute
            for andindex in range(len(nodeIndividual)):
                newval = ""
                if nodeIndividual[andindex] == 1:
                    # if a shadow and contributes, compute its value using its upstream nodes
                    if andNodeInvertList[andindex][0]:
                        newval = newval + "!"
                    newval = newval + self.nodeList[andNodes[andindex][0]]
                    for addnode in range(1, len(andNodes[andindex])):
                        newval = newval + " & "
                        if andNodeInvertList[andindex][addnode]:
                            newval = newval + " !"
                        newval = newval + self.nodeList[andNodes[andindex][addnode]]
                    orset.append(newval)
                # combine the shadow and nodes with or operations
            writenode = writenode + orset.pop()
            for val in orset:
                writenode = writenode + " | " + val
            # print(writenode)
            return writenode

    def writeModel_BoolNet(self, individual, model):
        # iterate over nodes to generate a BooleanNet representation for the entire model
        addString = ""
        # for i in range(0,len(self.nodeList)):
        for i in range(0, len(model.nodePositions)):
            addString = addString + model.writeNode_BoolNet(
                i,
                individual[model.individualParse[i] : model.individualParse[i + 1]],
                model,
            )
            addString = addString + "\n"
        return addString[:-1]

    class Params:
        def __init__(self):
            pass

        def simParams(
            self,
            mutModel=0.25,
            cells=1,
            samples=1,
            generations=10,
            popSize=24,
            mu=10,
            lambd=24,
            iters=100,
            genSteps=100,
            simSteps=100,
            crossoverProb=0.1,
            mutationProb=0.9,
            bitFlipProb=0.5,
        ):
            self.mutModel = mutModel
            self.cells = cells
            self.samples = samples
            self.generations = generations  # generations to run #100
            self.popSize = popSize  # size of population #24
            self.mu = mu  # individuals selected #24
            self.lambd = lambd  # children produced #24
            self.iters = iters  # number of simulations to try in asynchronous mode
            self.genSteps = genSteps  # steps to find steady state with fake data
            self.simSteps = (
                simSteps  # number of steps each individual is run when evaluating
            )
            self.crossoverProb = (
                crossoverProb  # prob of crossing over a particular parent
            )
            self.mutationProb = mutationProb  # prob of mutating a particular parent
            self.bitFlipProb = bitFlipProb  # prob of flipping bits inside mutation


class singleCell(ruleMaker):

    """Class for single-cell experiments"""

    def __init__(self, dataName, sep):
        """Read in pre-processed data and binarize by k-means"""
        data = np.loadtxt(dataName, delimiter=sep, dtype="str")
        self.geneList, self.sampleList, self.expMat = (
            data[1:, 0],
            data[0, 1:],
            sparse.csr_matrix(data[1:, 1:].astype("float")),
        )
        self.geneList = list(self.geneList)
        self.sampleList = list(self.sampleList)
        # print(self.sampleList)
        print("Genelist: " + str(self.geneList))
        # print(data)
        del data
        self.expMat.eliminate_zeros()
        self.binMat = preprocessing.binarize(
            self.expMat, threshold=0.001, copy=False
        )  # if data >=0: binMat = 1

        self.maxNodes = 15000  # 15000#
        self.maxSamples = 10000  # 60000#

        # self.binMat=self.expMat
        # super().__init__(self)
        # print(self.binMat.toarray())
        # populate binMat to a predefined size of 15000 genes and 10000 cells
        # print(self.binMat.shape)
        self.binMat.resize((self.maxNodes, self.maxSamples))
        # print(self.binMat)
        self.pathwayGraphs = {}
        # print(self.binMat.shape)
        # print(np.array(self.binMat))

    def addSubpop(self, subpopFile, sep):
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

    def filterData(self, threshold):
        """Filter data based on CV cutoffs"""
        self.cvGenes = []
        if threshold is not None:
            for i in range(0, self.expMat.get_shape()[0]):
                rowData = list(self.expMat.getrow(i).todense())
                if np.std(rowData) / np.mean(rowData) >= threshold:
                    self.cvGenes.append(self.geneList[i])
        else:
            self.cvGenes = copy.deepcopy(self.geneList)

    def find_pathways_from_gmt_LEGACY(self, gmtFile, minOverlap=4):
        """LEGACY CODE: find list of pathways with at least four genes found in data"""
        if not hasattr(self, cvGenes):
            print("You have not filtered the genes; are you sure you want to proceed?")
        self.pathwayGraphs = []
        if hasattr(self, "cvGenes"):
            pathwayGenes = set(self.cvGenes)
        else:
            if not hasattr(self, "cvGenes"):
                pathwayGenes = set(self.geneList)
        self.pathways = {}
        GMTdict = read_gmt(gmtFile)
        for key in GMTdict.keys():
            if (
                len(pathwayGenes.intersection(GMTdict[key])) >= minOverlap
            ):  # ensure there are at least minOverlap nodes in both pathway and detected genes
                print(key)
                self.pathways[key] = pathwayGenes.intersection(GMTdict[key])
                print(len(self.pathways[key]))
                # Get the KEGG network
                sif = makeSIF(key)
                pathGraph = sif_to_digraph(sif)
                self.pathwayGraphs.append(pathGraph)

    def find_pathways(
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

    def add_wikipathways(self, minOverlap=20, writeGraphml=True, removeSelfEdges=False):
        if hasattr(self, "cvGenes"):
            pathwayGenes = set(self.cvGenes)
        elif not hasattr(self, "cvGenes"):
            print("You have not filtered genes by any criterion.")
            pathwayGenes = set(self.geneList)

        for x in list(glob.glob("*[0-9].graphml")):
            G = nx.read_graphml(x)
            nodes = set(G.nodes())
            test = len(nodes.intersection(pathwayGenes))

            if test >= minOverlap:
                print("Pathway: ", x, " Overlap: ", test, " Edges: ", len(G.edges()))
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

    def makemetaNetwork(self):
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

    # make empty list representing no knockouts or knockins
    def setupEmptyKOKI(self):
        self.knockoutLists = [0] * len(self.nodePositions)
        self.knockinLists = [0] * len(self.nodePositions)
        # for q in range(len(self.sampleList)):
        # for q in range(len(self.nodePositions)):
        # self.knockoutLists.append([])
        ##self.knockoutLists[q]=0
        # self.knockinLists.append([])
        ##self.knockinLists[q]=0

    # calculate importance scores
    def calcImportance(
        self, individual, params, model, sss, knockoutLists, knockinLists, scSyncBoolC
    ):
        importanceScores = []  # holder for impact scores

        for node in self.nodePositions:
            print(node)
            """
            boolValues1 = self.evaluateByNode(individual, self.params.cells, model, self.nodePositions, self.params, node, [], scSyncBoolC) #knocked out
            boolValues2 = self.evaluateByNode(individual, self.params.cells, model, self.nodePositions, self.params, [], node, scSyncBoolC) #knocked in
            print(boolValues1)
            print(boolValues2)
            SSE=[abs(KO - KI) for KO, KI in zip(boolValues1, boolValues2)]
            print(SSE)
            SSE=sum(SSE)
            print(SSE)
            importanceScores.append(SSE)
            """
            # NP(self, individual, model, cells, params, KOs, KIs, scSyncBoolC)
            KO_temp = copy.deepcopy(self.knockoutLists)
            KO_temp[node] = 1
            print(KO_temp)
            boolValues1 = self.NP(
                individual,
                model,
                self.params.cells,
                self.params,
                KO_temp,
                knockinLists,
                scSyncBoolC,
            )  # knocked out
            KI_temp = copy.deepcopy(self.knockinLists)
            KI_temp[node] = 1
            print(KI_temp)
            boolValues2 = self.NP(
                individual,
                model,
                self.params.cells,
                self.params,
                knockoutLists,
                KI_temp,
                scSyncBoolC,
            )  # knocked in
            importanceScores.append(
                self.calculateImportanceScoreFromSimData(boolValues1, boolValues2)
            )
        print(importanceScores)

        """
        # loop over all nodes
        #for node in range(len(model.nodeList)):
        #for node in self.nodePositions:
        for node in range(0, len(self.nodePositions)):
            #SSEs=[] # add SSE across samples
            #nodeValues=[sss[j][model.nodeList[node]]for j in range(0,len(sss))]
            #nodeValues=self.binMat[node, range(0, len(sss))]
            #nodeValues=list(nodeValues.toarray()[0])
            #print(nodeValues)
            #return tuple([np.sum([(boolValues[j][i1]-self.binMat[i2, j])**2 for j in range(0,len(sss))]) for i2, i1 in zip(self.nodePositions, range(0, len(model.nodeList)))]) #fixed now
            #for j in range(0,len(sss)): # knock each node out and in to observe differences
            #for j in range(0, len(self.nodeList)):
            for j in range(0, len(self.nodePositions)):
                # print(j)
                #ss=sss[j]
                #initValues=nodeValues[j] #list(model.initValueList[j])
                # print(len(initValues))
                # print("#####")
                # print(initValues)
                # print("#####")
                knockerOuter=list(knockoutLists[j])
                # print(len(knockerOuter))
                # print("#####")
                knockerOuter.append(node) # knock out node
                # print(knockerOuter)
                # print("#####")
                # print(knockinLists[j])
                #boolValues1=self.NP(individual, model, self.params.cells, initValues, self.params, knockerOuter, knockinLists[j], scSyncBoolC) #def NP(self, individual, model, cells, sampleProbs, params, KOs, KIs, syncBoolC)
                print(knockerOuter)
                print(knockinLists[j])

                boolValues1 = self.evaluateByNode(individual, self.params.cells, model, self.nodePositions, self.params, knockerOuter, knockinLists[j], scSyncBoolC)

                knockerInner=list(knockinLists[j])
                knockerInner.append(node) # knock in node

                print(knockerInner)
                print(knockoutLists[j])
                boolValues2 = self.evaluateByNode(individual, self.params.cells, model, self.nodePositions, self.params, knockoutLists[j], knockerInner, scSyncBoolC)

                print(boolValues1)
                print(boolValues2)

                SSE=[abs(KO - KI) for KO, KI in zip(boolValues1, boolValues2)]
                print(SSE)
                SSE=sum(SSE)
                # find difference between knockout and knockin
                #SSE=0
                #for i in range(0, len(model.nodeList)):
                #for i in range(0, len(self.nodePositions)):
                #    SSE+=(boolValues1[i]-boolValues2[i])**2
                #SSEs.append(SSE)
            #importanceScores.append(sum(SSEs))
            #importanceScores.append(SSE)
        #print(importanceScores)
        """
        return importanceScores

    # calculate z score for a given pathway
    def scorePathway(self, RAs, pathImportances, pathname):
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
        plt.title(getPathwayName(str(pathname)))
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

    def run_pathway_analysis_from_outs(self, contrast, conditions, delimiter):
        self.binMat2 = self.binMat.A
        self.expMat2 = self.expMat.A
        pvalDict = {}
        overallUpreg = {}
        # open the set of differences to be considered
        contrastList = []
        for row in csv.reader(open(contrast), delimiter=delimiter):
            contrastList.append(row)
        # contrasts=contrasts[0] for testing purposes, assume that we have just one contrast (true in our AS+/AS- case study)
        conditions = pd.read_csv(conditions, sep=delimiter)
        print(conditions)
        # self.add_wikipathways(minOverlap=1)

        for contrasts in contrastList:
            condition1 = conditions[[str(contrasts[0])]]
            condition2 = conditions[[str(contrasts[1])]]
            # cells_condition1=(conditions[conditions[str(contrasts[0])] == 1])['Samples'].tolist()
            cells_condition1 = conditions.loc[
                conditions[str(contrasts[0])] == 1, list(conditions.columns)[0]
            ]
            cells_condition1 = list(cells_condition1)
            print(len(cells_condition1))
            # cells_condition2=(conditions[conditions[str(contrasts[1])] == 1])['Samples'].tolist()
            cells_condition2 = conditions.loc[
                conditions[str(contrasts[1])] == 1, list(conditions.columns)[0]
            ]
            cells_condition2 = list(cells_condition2)
            print(len(cells_condition1))
            index_condition1 = [
                self.sampleList.index(i)
                for i in set(cells_condition1).intersection(set(self.sampleList))
            ]
            index_condition2 = [
                self.sampleList.index(i)
                for i in set(cells_condition2).intersection(set(self.sampleList))
            ]

            # print(cells_condition1)
            # print(len(set(cells_condition1).intersection(set(self.sampleList))))
            # print(len(set(cells_condition2).intersection(set(self.sampleList))))

            # make RA - in the case of single cell experiments, find the proportion of cells in which the gene is expressed
            for pathname in list(self.pathwayGraphs.keys()):
                print(pathname)

                if os.path.exists(pathname[:-8] + "_IS.graphml") or os.path.exists(
                    pathname[:-8] + ".graphml_processed.graphml_importanceScores.csv"
                ):
                    """
                    if os.path.exists(pathname[:-8]+"_IS.graphml"):
                        print(str(pathname[:-8]+"_IS.graphml"))

                        paNetTemp = nx.read_graphml(pathname[:-8]+"_IS.graphml")
                        # get nodeScores
                        nodeScores = dict(paNetTemp.nodes(data='importanceScore', default=np.nan))

                        # Calculate RA
                        RA={}
                        for node in list(paNetTemp.nodes()):
                            node_index=self.geneList.index(node)
                            expression_condition1=np.mean(self.binMat2[node_index, index_condition1].tolist())
                            expression_condition2=np.mean(self.binMat2[node_index, index_condition2].tolist())
                            RA[node]=abs(expression_condition1-expression_condition2)

                        # add RA as attribute to graph
                        nx.set_node_attributes(paNetTemp, values=RA, name='relativeAbundance')

                        z_scores=[]
                        # iterate over comparisons for each pathway and calculate z score
                        zscore, impactScore = self.scorePathway(RA, nodeScores)
                        z_scores.append(zscore)
                        pvals = scipy.stats.norm.sf(z_scores) # calculate p value
                        #print(pvals)
                        pvalDict[str(pathname)]=pvals

                        # add impact score as attribute to graph
                        nx.set_node_attributes(paNetTemp, values=impactScore, name='impactScore')

                        # write out graph with additions
                        nx.write_graphml_lxml(paNetTemp, pathname[:-8]+"_IS_"+"_vs_".join(contrasts)+".graphml")
                    else:
                    """
                    if os.path.exists(
                        pathname[:-8]
                        + ".graphml_processed.graphml_importanceScores.csv"
                    ):
                        print(
                            pathname[:-8]
                            + ".graphml_processed.graphml_importanceScores.csv"
                        )
                        paNetTemp = nx.read_graphml(
                            pathname[:-8] + ".graphml_processed.graphml"
                        )
                        # get nodeScores
                        nodeScoresDF = pd.read_csv(
                            pathname[:-8]
                            + ".graphml_processed.graphml_importanceScores.csv"
                        )
                        nodeScoresDF.index = list(nodeScoresDF.Node)
                        if "Strat3_IS" in nodeScoresDF.columns:
                            # add impact score as attribute to graph
                            nodeScores = dict(
                                zip(
                                    list(nodeScoresDF.index),
                                    list(nodeScoresDF.loc[:, "Strat3_IS"]),
                                )
                            )
                            nx.set_node_attributes(
                                paNetTemp, values=nodeScores, name="Strat3_IS"
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
                            nx.set_node_attributes(
                                paNetTemp, values=obsERS, name="Max ERS"
                            )
                            # write out IS graph with additions
                            # nx.write_graphml_lxml(paNetrTemp, pathname[:-8]+"_IS.graphml")
                            # Calculate RA
                            RA = {}
                            upreg_condition1 = {}
                            expression_condition1 = {}
                            expression_condition2 = {}
                            binexpression_condition1 = {}
                            binexpression_condition2 = {}
                            nodeScoresDF[str(contrasts[0])] = np.nan
                            nodeScoresDF[str(contrasts[1])] = np.nan
                            nodeScoresDF["BIN_" + str(contrasts[0])] = np.nan
                            nodeScoresDF["BIN_" + str(contrasts[1])] = np.nan
                            for node in list(nodeScoresDF.index):
                                node_index = self.geneList.index(node)
                                binexpression_condition1[node] = np.mean(
                                    self.binMat2[node_index, index_condition1].tolist()
                                )
                                binexpression_condition2[node] = np.mean(
                                    self.binMat2[node_index, index_condition2].tolist()
                                )
                                expression_condition1[node] = np.mean(
                                    self.expMat2[node_index, index_condition1].tolist()
                                )
                                expression_condition2[node] = np.mean(
                                    self.expMat2[node_index, index_condition2].tolist()
                                )
                                print(
                                    [
                                        expression_condition1[node],
                                        expression_condition2[node],
                                    ]
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
                            nodeScoresDF["BIN_" + str(contrasts[0])] = nodeScoresDF[
                                "Node"
                            ].map(binexpression_condition1)
                            nodeScoresDF["BIN_" + str(contrasts[1])] = nodeScoresDF[
                                "Node"
                            ].map(binexpression_condition2)
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
                            zscore, impactScore, CVdict = self.scorePathway(
                                modRA, nodeScores, pathname[:-8]
                            )
                            z_scores.append(zscore)
                            pvals = scipy.stats.norm.sf(z_scores)  # calculate p value
                            # print(pvals)
                            pvalDict[str(pathname)] = [
                                pvals,
                                str(CVdict),
                                str(zscore),
                                str(sum(impactScore.values())),
                                str(mean(impactScore.values())),
                            ]

                            # add impact score as attribute to graph
                            # nx.set_node_attributes(paNetTemp, values=impactScore, name='impactScore')

                            # write out graph with additions if pval < 0.05
                            if pvalDict[str(pathname)][0] < 0.05:
                                nx.write_graphml_lxml(
                                    paNetTemp,
                                    pathname[:-8]
                                    + "_IS_"
                                    + "_vs_".join(contrasts)
                                    + ".graphml",
                                )
                            nodeScoresDF.to_csv(
                                pathname[:-8]
                                + ".graphml_processed.graphml_importanceScores.csv",
                                index=False,
                            )
                            if (
                                sum(
                                    nodeScoresDF[
                                        str("Upregulated_in_" + str(contrasts[0]))
                                    ]
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
                # print(key, str(getPathwayName(key)), str(value[0]))
                # fh.write(','.join([key, str(getPathwayName(key[:8])), str(value[0]),"_vs_".join(contrasts), str(overallUpreg[key[:8]])]))
                # fh.write(','.join([key, str(getPathwayName(key[:8])), str(value[0]), "_vs_".join(contrasts)])) #, str(overallUpreg[key])]))
                pvalDict[key] = [
                    key,
                    str(getPathwayName(key[:8])),
                    str(value[0][0]),
                    "_vs_".join(contrasts),
                    str(overallUpreg[key]),
                    str(value[1]),
                    str(value[2]),
                    str(value[3]),
                    str(value[4]),
                ]
                # fh.write("\n")
            # fh.close()
            pvalDF = pd.DataFrame.from_dict(pvalDict, orient="index")
            pvalDF.columns = [
                "Pathway ID",
                "Pathway Name",
                "P value",
                "Contrast",
                str("Upregulated_in_" + str(contrasts[0])),
                "CVdict",
                "zscore",
                "impactScore",
                "meanNodewiseImpactScore",
            ]
            pvalDF.to_csv("pvalues_" + "_vs_".join(contrasts) + ".csv", index=False)
        del self.binMat2

    def scoreNodes(self, graph):
        """Skeleton for node scoring"""
        # Notes: please disentangle ruleMaker and singleCell classes, currently it looks like ruleMaker methods are being indiscriminately called on singleCell
        # model={}
        # Task 1: Retrieve data for the pathway being scored
        # a. get nodes
        # Read in a graphml file. Assume that the wrapper for this function is going through each graphml file in the directory.
        net = nx.read_graphml(graph)
        # nodes=list(net)
        # model["Network"]=net
        # b. subset binMat
        netGenes = [self.geneList.index(gene) for gene in list(net)]
        # c. store data in accessible form
        # data=scObj.binMat[netGenes, :]
        # print(data.shape)
        # print(data.todense())
        # d. provide for destruction of data structure
        # del data
        # Task 2: Narrow down predecessors of each node
        # Go to function getPreds
        # self.getPreds(geneList, net)
        # Task 3: Infer rules for each subpopulation individually
        # model["initValueList"]=self.__genInitValueList(graph=net)
        # print(model.keys())
        # read in C function to run simulations
        # updateBooler=ctypes.cdll.LoadLibrary('./simulator.so')
        # set up parameters of run, model
        # ruleMaker is the parent class of singleCell; super().method(args) calls methods from the super class. For instance: super().__init__(net) calls the super constructor from ruleMaker to self, ie, on the singleCell object. Similarly for simParams and updateCpointers, both of which are defined in the ruleMaker class
        # super().simParams()
        updateBooler = ctypes.cdll.LoadLibrary("./simulator.so")
        scSyncBoolC = updateBooler.scSyncBool
        super().__init__(net)
        super().updateCpointers()
        self.setupEmptyKOKI()
        # self.__genInitValueList(net)
        print("Node positions:\t" + str(self.nodePositions))
        # print(self.nodePositions)
        # find rules by doing GA then local search
        # GA
        population, logbook = self.eaMuPlusLambdaAdaptive(scSyncBoolC, graph)

        out1, out2, model = self.findPopBest(population)

        print(self.nodePositions)
        print(model)
        print(out1)
        print(out2)

        with open(graph + "_rules_GA.txt", "w") as text_file:
            text_file.write(model.writeModel(out2, model))

        pickle.dump(out2, open(graph + "_out2.pickle", "wb"))

        # model.setupEmptyKOKI()
        # Local search
        outputs = [
            self.checkNodePossibilities(
                node,
                out2,
                model.sampleList,
                model.params.cells,
                model,
                model.params,
                model.knockoutLists,
                model.knockinLists,
                scSyncBoolC,
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
        # write rules
        pickle.dump(equivs, open(graph + "_equivs1.pickle", "wb"))
        # write local errors
        pickle.dump(localErrors, open(graph + "_localErrors1.pickle", "wb"))

        bruteout2 = []
        for i in range(len(equivs)):
            bruteout2.extend(equivs[i][randint(0, len(equivs[i]) - 1)])

        with open(graph + "_rules_LS.txt", "w") as text_file:
            text_file.write(model.writeModel(bruteout2, model))

        # pickle.dump(model, open(graph+"_model.pickle", "wb"))
        pickle.dump(bruteout2, open(graph + "_bruteout2.pickle", "wb"))

        # calculate importance scores and output
        print(self.knockinLists)
        print(self.knockoutLists)
        scores1 = self.calcImportance(
            bruteout2,
            model.params,
            model,
            model.sampleList,
            model.knockoutLists,
            model.knockinLists,
            scSyncBoolC,
        )
        pickle.dump(scores1, open(graph + "_scores1.pickle", "wb"))
        # print out importance scores
        fh = open(graph + "_importance_scores.txt", "w+")
        importanceVals = {}
        for node in range(0, len(list(self.nodeList))):
            node1 = self.nodeList[node]
            importanceVals[node1] = scores1[node]
            print(node1, scores1[node])
            fh.write(",".join([str(node1), str(scores1[node])]))
            fh.write("\n")
        fh.close()
        nx.set_node_attributes(net, importanceVals, "nodeScore")
        nx.write_graphml(net, graph + "_scored.graphml")

    def inherit(
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


def read_gmt(filename):
    """read in file with pathway genes and names"""
    gmt_dict = {}
    inputfile = open(filename, "r")
    lines = inputfile.readlines()
    for line in lines:
        newline = line.split("    ")
        gmt_dict[newline[0]] = set(newline[2:])
    return gmt_dict


def xstr(s):
    """Handle NoneType in strings"""
    if s is None:
        return "Unknown"
    return str(s)


def hms_string(sec_elapsed):
    # write output as hours:minutes:seconds
    h = int(sec_elapsed / (60 * 60))
    m = int((sec_elapsed % (60 * 60)) / 60)
    s = sec_elapsed % 60.0
    return "{}:{:>02}:{:>05.2f}".format(h, m, s)


def supplementBONITA(pathwayGraphs):

    # Get number of test networks:
    print("Number of test networks: ", len(pathwayGraphs.keys()))
    # print("Pathway Name,", "Length of shortest longest path,","Number of nodes in pathway,","Source nodes,","Number of source nodes,","Out-degree of source nodes,", "Mean out-degree of source nodes")
    fh = open("testSimParams.csv", "w+")
    fh.write(
        "".join(
            [
                "Pathway ID,",
                "Pathway Name,",
                "Length of shortest longest path,",
                "Number of nodes in pathway,",
                "Number of source nodes,",
                "Mean out-degree of source nodes",
                "\n",
            ]
        )
    )
    # For all pathway graphs
    for name, path in pathwayGraphs.items():
        # Get pathway name
        pathName = getPathwayName(name)
        print(pathName)
        # Find number of nodes
        numberNodes = nx.number_of_nodes(path)
        # Find source nodes
        sourceNodes = []
        for x in path.in_degree():
            if x[1] == 0:
                sourceNodes.append(x[0])
        # Find out-degree of source nodes
        sourceOutDeg = []
        for x in sourceNodes:
            sourceOutDeg.append(path.out_degree(x))
        # Find mean out-degree of source nodes
        if len(sourceOutDeg) >= 1:
            sourceOutDeg_mean = mean(sourceOutDeg)
        else:
            sourceOutDeg_mean = "Not defined"
        # Find longest shortest path
        tempMax = -1
        temp = dict(nx.all_pairs_shortest_path_length(path, cutoff=None))
        for x, y in temp.items():
            if max(y.values()) >= tempMax:
                tempMax = max(y.values())
        # print("Pathway: ", name, "Length of shortest longest path: ", tempMax, "Number of nodes in pathway: ", numberNodes, "Source nodes: ", sourceNodes, "Number of source nodes: ", len(sourceNodes), "Out-degree of source nodes: ", sourceOutDeg, "Mean out-degree of source nodes: ", sourceOutDeg_mean)
        # print(name, ",", tempMax, ",", numberNodes, ",", sourceNodes, ",", len(sourceNodes),",", sourceOutDeg, ",", sourceOutDeg_mean)
        fh.write(
            "".join(
                [
                    str(name),
                    ",",
                    str(pathName),
                    ",",
                    str(tempMax),
                    ",",
                    str(numberNodes),
                    ",",
                    str(len(sourceNodes)),
                    ",",
                    str(sourceOutDeg_mean),
                    "\n",
                ]
            )
        )
    fh.close()


def pipeline2(writeGraphml=False, organism="hsa"):
    pathwayGraphs = {}
    k = KEGG()
    k.organism = organism
    for x in list(k.pathwayIds):
        print(x)
        x = x.replace("path:", "")
        # x=x.replace("hsa","ko")
        sif = makeSIF(pathway=x, keggObject=k)
        print(len(sif))
        if len(sif) < 1:
            i = 0
        else:
            nodes = set(sif[0:, 0]).union(set(sif[0:, 2]))
            print("Pathway: ", x, " Edges: ", len(sif))
            G = sif_to_digraph2(sif=sif)
            pathwayGraphs[x] = G
            if writeGraphml:
                nx.write_graphml(G, x + ".graphml", infer_numeric_types=True)

    supplementBONITA(pathwayGraphs)
    # import pandas
    # import seaborn as sns
    # resFile=pandas.read_csv("testSimParams.csv")
    # print(resFile)
    # sns.set_palette(sns.color_palette("Greys_r", 6))
    # sns.set_context(context='paper', font_scale=0.8)
    # sns.set_style("white")
    # b_width = 1  # chose an arbitrary value here
    # my_bins = np.arange(min(temp_rnAllNodes), max(temp_rnAllNodes) + b_width, b_width)
    # fig2=sns.distplot(temp_rnAllNodes, color='grey', kde=False, bins=my_bins)
    # plt.xlim(1,7)
    # fig2.set(xlabel='Rule Number')
    # fig2=fig2.get_figure()
    # fig2.savefig("Overall_rnAllnodes_dist.svg")
    # plt.close()


def pipeline(dataName="", sep=",", subpopFile=""):

    print("##########")
    print("Create singleCell object from matrix of normalized read counts.")
    scTest = singleCell(dataName=dataName, sep=sep)
    print("##########")
    print("Dimensions of test dataset: ", scTest.binMat.shape)
    print("##########")
    # #print("Read in file with subpopulation information and add information to singleCell object.")
    # #scTest.addSubpop(subpopFile="subpopulations.csv", sep=",")
    # print("##########")
    print("Filter data using CV >=0")
    print("##########")
    scTest.filterData(threshold=None)
    print("Dimensions of test dataset: ", len(scTest.cvGenes))
    print("##########")

    print("Find pathways with at least 25 genes in dataset")

    scTest.add_wikipathways(minOverlap=1)

    scTest.find_pathways(
        minOverlap=10,
        predefinedList=[
            "hsa04620",  # Toll-like receptor signaling pathway
            "hsa04650",  # Natural killer cell mediated cytotoxicity
            "hsa04640",  # Hematopoietic cell lineage
            "hsa04630",  # JAK-STAT signaling pathway
            "hsa04350",  # TGF-beta signaling pathway
            "hsa04110",  # Cell cycle
            "hsa04610",  # Complement and coagulation cascades
            "hsa04611",  # Platelet activation
            "hsa04624",  # Toll and Imd signaling pathway
            "hsa04621",  # NOD-like receptor signaling pathway
            "hsa04622",  # RIG-I-like receptor signaling pathway
            "hsa04623",  # Cytosolic DNA-sensing pathway
            "hsa04625",  # C-type lectin receptor signaling pathway
            "hsa04612",  # Antigen processing and presentation
            "hsa04660",  # T cell receptor signaling pathway
            "hsa04658",  # Th1 and Th2 cell differentiation
            "hsa04659",  # Th17 cell differentiation
            "hsa04657",  # IL-17 signaling pathway
            "hsa04662",  # B cell receptor signaling pathway
            "hsa04664",  # Fc epsilon RI signaling pathway
            "hsa04666",  # Fc gamma R-mediated phagocytosis
            "hsa04670",  # Leukocyte transendothelial migration
            "hsa04672",  # Intestinal immune network for IgA production
            "hsa04062",  # Chemokine signaling pathway
        ],
    )

    # cut down data size
    nodeIndices = []
    for graph in scTest.pathwayGraphs.keys():
        graph = scTest.pathwayGraphs[graph]
        # collect list of genes from pathways found above
        print(graph.nodes())
        # find indices of these genes in self.geneList
        for node in graph.nodes():
            nodeIndices.append(scTest.geneList.index(node))

    nodeIndices = set(nodeIndices)
    nodeIndices = list(nodeIndices)
    print(nodeIndices)
    print(len(nodeIndices))

    # retain only rows of those genes from binMat
    scTest.expMat = scTest.expMat[nodeIndices, :]
    print(scTest.expMat.shape)
    scTest.binMat = preprocessing.binarize(
        scTest.expMat, threshold=0.001, copy=False
    )  # if data >=0: binMat = 1
    scTest.maxNodes = 15000  # 25000#200
    scTest.maxSamples = 10000  # 60000#10000
    print(scTest.binMat.shape)
    print(scTest.binMat.A[0:10, 0:10])
    print((scTest.maxNodes, scTest.maxSamples))
    scTest.binMat.resize(
        (scTest.maxNodes, scTest.maxSamples)
    )  # restrained to 1000 genes and 1000 samples for testing simulation code, so that I can have a 200 step simulation
    print(scTest.binMat.A[0:10, 0:10])
    print(scTest.binMat.shape)
    print("Old genelist:\t" + str(scTest.geneList))
    scTest.geneList = [scTest.geneList[node] for node in nodeIndices]
    print("New genelist:\t" + str(scTest.geneList))
    scTest.nodeList = scTest.geneList
    # scTest.nodeList.sort()
    scTest.nodePositions = [scTest.geneList.index(node) for node in scTest.nodeList]

    print("##########")
    # print("Prepare a combined network.")
    # scTest.makemetaNetwork()
    # nx.write_graphml(scTest.metaNetwork, dataName+"metanetwork.graphml", infer_numeric_types=True)
    pickle.dump(scTest, open(dataName + "scTest.pickle", "wb"))
    # scTest.supplementBONITA() # for bulk bonita supplement
    # plt.subplot(121)
    # nx.draw_shell(testG, with_labels=True, font_weight='bold')
    # plt.show()
    # plt.subplot(122)
    # nx.draw(testG, pos=nx.spring_layout(testG))
    # plt.show()

    # print("##########")
    # print("Find correlations between each node and its predecessors.")
    # corHist1=[]
    # for graph in glob.glob('*.graphml'):
    #     print("##########")
    #     print(graph)
    #     scTest.scoreNodes(graph)
    #     corHist1.append(scTest.corHist)

    # print(len(corHist1))
    # flattened_list = [y for x in corHist1 for y in x]
    # print(len(flattened_list))
    # corHistPlot=seaborn.distplot(flattened_list , bins=None, hist=True, kde=False, rug=True, fit=None, hist_kws=None, kde_kws=None, rug_kws=None, fit_kws=None, color="g", vertical=False, norm_hist=False, axlabel="Mutual information to predecessor nodes", label=None, ax=None)
    # fig=corHistPlot.get_figure()
    # fig.savefig(dataName+"corHistPlot.png")
    # print("##########")
    # print("Plotted histogram of mutual information between nodes and their predecessors - corHistPlot.png")
    # print("Mean of MI: ", mean(flattened_list))
    # print("Standard deviation of MI: ", stdev(flattened_list))
    # print("##########")
    # print("Set up simulation parameters.")

    scTest = pickle.load(open(dataName + "scTest.pickle", "rb"))
    runAllNets()
    # scTest.scoreNodes(graph="hsa05418.graphml")
    # scTest.scoreNodes(graph="allscData.txtmetanetwork.graphml")
    # print("Execution time: {}\n\n".format(hms_string(time() - start_time)))


def runAllNets():
    for net in glob.glob("*_processed.graphml"):
        name = net[:-8]
        shellHandle = open(name + "_scoreNodes.sh", "w+")
        # slurmCommands=str("#!/bin/sh\n#SBATCH --partition=standard\n#SBATCH -J "+name+"\n#SBATCH -o "+name+".out\n#SBATCH -t 72:00:00\n#SBATCH -n 1\n#SBATCH -c 1\n#SBATCH --mem=100G\nmodule load intelpython3/2019.3\npython3.6 scylight_spawn_pathwayAnalysis.py "+"False "+str(net))
        # slurmCommands=str("#!/bin/sh\n#SBATCH --partition=debug\n#SBATCH -J "+name+"\n#SBATCH -o "+name+".out\n#SBATCH -t 1:00:00\n#SBATCH -n 1\n#SBATCH -c 1\n#SBATCH --mem=50G\nmake\nmodule load intelpython3/2019.3\npython3.6 scylight_spawn_pathwayAnalysis.py "+"False "+str(net))
        # slurmCommands=str("#!/bin/sh\n#SBATCH --partition=preempt\n#SBATCH -J "+name+"\n#SBATCH -o "+name+".log\n#SBATCH -t 10:00:00\n#SBATCH -n 1\n#SBATCH -c 1\n#SBATCH --mem=120G\nmodule load intelpython3/2019.3\nmake\npython3.6 scylight_spawn_pathwayAnalysis.py "+"False "+str(net))
        slurmCommands = str(
            "#!/bin/sh\n#SBATCH --partition=standard\n#SBATCH -J "
            + name
            + "\n#SBATCH -o "
            + name
            + ".log\n#SBATCH -t 10:00:00\n#SBATCH -n 1\n#SBATCH -c 1\n#SBATCH --mem=120G\nmodule load intelpython3/2019.3\nmake\npython3.6 scylight_spawn_pathwayAnalysis.py "
            + "False "
            + str(net)
        )
        shellHandle.write(slurmCommands)
        shellHandle.close()
        shellCommand = name + "_scoreNodes.sh"
        print([shellCommand])
        p = subprocess.Popen(["sbatch", shellCommand])


# Test
# for subpop in glob.glob("*.bin"):
#    print(subpop)
#    pipeline(dataName=subpop, sep=",")

# pipeline(dataName="simstudy.txt", subpopFile="subpopulations.csv")
# pipeline(dataName="allscData.txt", subpopFile="subpopulations.csv")


if __name__ == "__main__":

    start_time = time.time()

    # read in arguments from shell scripts
    parser = argparse.ArgumentParser()
    parser.add_argument("fullPipeline")
    parser.add_argument("network")
    results = parser.parse_args()
    fullPipeline = results.fullPipeline
    net = results.network

    if fullPipeline == "True":
        dataFile = glob.glob("*.bin")[0]
        print(dataFile)
        pipeline(dataName=dataFile, subpopFile="subpopulations.csv", sep=",")
        # pipeline(dataName="allscData.bin", subpopFile="subpopulations.csv", sep=",")
    else:
        dataFile = glob.glob("*.bin")[0]
        print(dataFile)
        scTest = pickle.load(open(dataFile + "scTest.pickle", "rb"))
        scTest.scoreNodes(graph=str(net))

# Test pathway analysis
# scTest = pickle.load(open("scTest.pickle", "rb"))
# scTest.run_pathway_analysis_from_outs("contrasts.txt", "conditions.txt")

# Test subpopulation selection
# scTest = pickle.load(open("scTest.pickle", "rb"))
# scTest.addSubpop("subpop_info.txt", "\t")
# print(set(list(scTest.subpopInfo.values())))

# reverseSubpopDict={}
# for subpop in set(list(scTest.subpopInfo.values())):
#     reverseSubpopDict[subpop]=[]

# for cell in list(scTest.subpopInfo.keys()):
#     tempSubpop=scTest.subpopInfo[cell]
#     reverseSubpopDict[tempSubpop].append(cell)

# import os
# for key, value in reverseSubpopDict.items():
#     #key = subpop; value = list of cells
#     cell_index=[i for i, e in enumerate(scTest.sampleList) if e in value]
#     tempBinMat=scTest.binMat[:, cell_index] #.todense()
#     tempBinMat=pandas.DataFrame(tempBinMat.toarray())
#     tempBinMat.to_csv(path_or_buf=key+"_subpop.csv", sep='\t')
#     os.mkdir(key+"_subpop")
