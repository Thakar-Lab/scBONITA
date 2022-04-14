from deap import base, creator, gp, tools
from deap import algorithms as algo
import numpy as np
import networkx as nx
from sklearn import preprocessing
from scipy.stats.stats import spearmanr
import ctypes as ctypes
import itertools as itertool
import copy
import pickle
from random import random, randint, sample, choice
import math
from collections import defaultdict
from itertools import chain
from operator import attrgetter
import gc
import pandas as pd


class ruleMaker:
    def __makeToolBox(self, graph):
        """sets up GA toolbox from deap"""
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
        toolbox.register("genRandomBitString", self.__genBits)  # , model=self)
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
        toolbox.register("select", self.__selNSGA2)
        toolbox.register("similar", np.array_equal)

        # ADD TOOLBOX TO OBJECT
        self.toolbox = toolbox
        self.stats = stats

    def __init__(
        self,
        graph,
        removeSelfEdges=False,
        restrictIncomingEdges=True,
        maxIncomingEdges=3,
        groundTruth=False,
        graphName="",
    ):
        """Initialize a ruleMaker object for rule inference with scBONITA - RD"""
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
            self.geneList.index(node) for node in nodeList if node in self.geneList
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
        self.possibilityList = possibilityLister
        self.possibilityInverter = possibilityInverter
        self.nodeNum = len(nodeList)
        self.params = self.Params()
        self.params._Params__simParams()
        self.__makeToolBox(graph)
        self.ruleGraph = ruleGraph
        self.nodeDict = nodeDict  # identifies names of nodes with their index in the node list.. provide name, get index
        self.successorNums = succnum
        # nx.write_graphml(ruleGraph, graphName+"_ruleGraph.graphml")
        print("\nIndividual parse: " + str(self.individualParse))
        print("\nNodelist: " + str(self.nodeList))
        print("\nNode positions: " + str(self.nodePositions))
        print("\nPossibilityList: " + str(self.possibilityList))

    def __update_upstream(self, node, newUpstreams):
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

    def __updateCpointers(self):
        """set up C pointers with correct lengths to pass to simulation software in C"""
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
        #self.andNodeInvert = np.array(tempandinverter, dtype=np.intc, order="C")
        self.andNodeInvert = np.array(tempandinverter, dtype=object, order="C")
        #self.andNodes = np.array(tempandnoder, dtype=np.intc, order="C")
        self.andNodes = np.array(tempandnoder, dtype=object, order="C")

    def __genRandBits(self):
        """generates a random bitstring"""
        arr = np.random.randint(2, size=(int(self.size),))
        return list(arr)

    def __findEnd(self, node):
        if node == len(self.nodeList) - 1:
            end = self.size
        else:
            end = self.individualParse[node + 1]
        return end

    def __cxTwoPointNode(self, ind1, ind2):
        """Executes a two-point crossover on the input :term:`sequence`
        individuals. The two individuals are modified in place and both keep
        their original length.
        :returns: A tuple of two individuals.
        This function uses the :func:`~random.randint` function from the Python
        base :mod:`random` module.

        Modified from deap to cross over between rules = needed to account for bistring only being one of two components of individual
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
        ind1[0]._ruleMaker__updateCpointers()
        ind2[0]._ruleMaker__updateCpointers()
        return ind1, ind2

    def __findPopBest(self, population):
        """finds the lowest error individual in a population"""
        saveVal = -1
        minny = float("Inf")
        for i in range(len(population)):
            if np.sum(population[i].fitness.values) < minny:
                minny = np.sum(population[i].fitness.values)
                saveVal = i
        ultimate = population[saveVal]
        minvals = population[saveVal].fitness.values
        return minvals, ultimate[1], ultimate[0]

    def __NP(self, individual, model, cells, params, KOs, KIs, scSyncBoolC):
        """NP simulation code for synchronous simulation"""
        cellArray = []

        # set up knockin and knockout lists
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
        individualParse = np.array(model.individualParse, dtype=np.intc, order="C")
        andLenList = np.array(model.andLenList, dtype=np.intc, order="C")
        nodePositions1 = model.nodePositions
        nodePositionsC = np.array(nodePositions1, dtype=np.intc, order="C")
        simSteps = self.params.simSteps
        lenSamples1 = len(model.sampleList)
        binMatC1 = self.binMat.toarray(order="C")
        binMatC3 = np.transpose(
            np.array(copy.deepcopy(binMatC1), order="C", dtype=np.intc)
        )
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

        return vals

    def __varOrAdaptive(
        self, population, toolbox, lambda_, cxpb, mutpb, genfrac, mutModel
    ):
        """generates list of offspring to be compared... decides to do crossover or mutation"""
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
                ind1, ind2 = self.__cxTwoPointNode(ind1, ind2)
                del ind1.fitness.values
                offspring.append(ind1)
            elif op_choice < cxpb + mutpb:  # Apply mutation
                ind = toolbox.clone(choice(population))
                (ind,) = self.__mutFlipBitAdapt(ind, genfrac, mutModel)
                del ind.fitness.values
                offspring.append(ind)
            else:  # shouldn't happen... clone existing individual
                offspring.append(choice(population))
        return offspring

    def __selectMutNode(self, errors):
        """select node to mutate"""
        normerrors = [
            1.0 * error / np.sum(errors) for error in errors
        ]  # normalize errors to get a probability that the node  is modified
        probs = np.cumsum(normerrors)
        randy = random()  # randomly select a node to mutate
        return next(i for i in range(len(probs)) if probs[i] > randy)

    def __mutFlipBitAdapt(self, indyIn, genfrac, mutModel):
        """mutation algorithm"""
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
            focusNode = self.__selectMutNode(pseudoerrors)
        else:
            # if errors are relatively high, focus on nodes that fit the worst and have highest in-degree
            # calculate probabilities for mutating each node
            for i in range(len(errors)):
                temper = model.successorNums[i]
                if temper == 0:
                    errors[i] = errors[i] * len(model.possibilityList[i])
                else:
                    errors[i] = errors[i] * len(model.possibilityList[i]) * temper
            focusNode = self.__selectMutNode(errors)
        # perform mutation
        if model.andLenList[focusNode] > 1:
            # find ends of the node of interest in the individual
            start = model.individualParse[focusNode]
            end = model._ruleMaker__findEnd(focusNode)
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
                        # print(addNoder)
                    else:
                        recalc = np.cumsum([1.0 * rval / tempsum for rval in rvals])
                        # print(recalc)
                        addNoder = next(
                            i for i in range(len(recalc)) if recalc[i] > randy
                        )
                        # print(addNoder)
                    temppermup.append(upstreamAdders.pop(addNoder))
                    # print(rvals)
                    rvals.pop(addNoder)
                model._ruleMaker__update_upstream(focusNode, temppermup)
                model._ruleMaker__updateCpointers()
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

    def __genBits(self):
        # generate random bitlist
        startInd = list(self.__genRandBits())
        counter = 0
        # make sure bitlist isn't zero
        while np.sum(startInd) == 0 and counter < float("Inf"):
            startInd = list(self.__genRandBits())
            counter += 1
        # go through nodes and make sure that there are 1-5 ones in the random list
        for node in range(0, len(self.nodeList)):
            end = self.__findEnd(node)
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

    def __sortNondominatedAdapt(self, individuals, k, first_front_only=False):
        """
        Taken from deap and modified slightly to make pareto sorting less strict
        Sort the first *k* *individuals* into different nondomination levels
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
                if self.__dominated(fit_i, fit_j):
                    dominating_fits[fit_j] += 1
                    dominated_fits[fit_i].append(fit_j)
                elif self.__dominated(fit_j, fit_i):
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

    def __dominated(self, ind1, ind2):
        """TTaken from deap and modified slightly to make pareto sorting less strict.
        Return true if each objective of *self* is not strictly worse than
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

    def __assignCrowdingDist(self, individuals):
        """taken from deap. Assign a crowding distance to each individual's fitness. The
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

    def __selNSGA2(self, individuals, k):

        """Calculate  fitness for an individual. NSGA2 selection taken from deap
        Apply NSGA-II selection operator on the *individuals*. Usually, the
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
        pareto_fronts = self.__sortNondominatedAdapt(individuals, k)
        for front in pareto_fronts:
            self.__assignCrowdingDist(front)

        chosen = list(chain(*pareto_fronts[:-1]))
        k = k - len(chosen)
        if k > 0:
            sorted_front = sorted(
                pareto_fronts[-1], key=attrgetter("fitness.crowding_dist"), reverse=True
            )
            chosen.extend(sorted_front[:k])

        return chosen

    def __bitList(self, n, x):
        templist = [1 if digit == "1" else 0 for digit in bin(n)[::-1]]
        while len(templist) < x:
            templist.append(0)
        while (len(templist)) > x:
            templist.pop()
        return templist

    def writeModel(self, individual, model):
        """iterate over nodes to generate a BooleanNet representation for the entire model"""
        addString = ""
        for i in range(0, len(model.nodePositions)):
            addString = addString + model._ruleMaker__writeNode(
                i,
                individual[model.individualParse[i] : model.individualParse[i + 1]],
                model,
            )
            addString = addString + "\n"
        return addString[:-1]

    def __findInEdges(self, model, node):
        """find the incoming edges to each 'and' connection for a given node"""
        inEdges = []
        for lister in model.andNodeList[node]:
            tempTup = tuple(lister)
            inEdges.append(set(tempTup))
        return inEdges

    def __simplifyRule(self, rule, inEdges):
        """find the simplest form of a rule"""
        for i in range(len(rule)):
            if rule[i] == 1:
                for j in range(len(rule)):
                    if rule[j] == 1 and not i == j:
                        if inEdges[i].issubset(inEdges[j]):
                            rule[j] = 0
        return rule

    def __writeNode(self, currentNode, nodeIndividual, model):
        """write out evaluation instructions in BooleanNet format. This follows the exact same code as updateNode (for switch=0), but writes a string instead of actually updating the values of the nodes"""
        andNodes = model.andNodeList[
            currentNode
        ]  # find the list of shadow and nodes we must compute before computing value of current nodes
        andNodeInvertList = model.andNodeInvertList[
            currentNode
        ]  # find list of lists of whether input nodes need to be inverted (corresponds to inputOrder)
        writenode = (
            "" + model.nodeList[currentNode] + " *= "
        )  # set up the initial string to use to write node
        inEdges = self.__findInEdges(model, currentNode)
        nodeIndividual = self.__simplifyRule(nodeIndividual, inEdges)

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

    def __writeNode_BoolNet(self, currentNode, nodeIndividual, model):
        """write out evaluation instructions in BoolNet format.
        This follows the exact same code as updateNode (for switch=0), but writes a string instead of actually updating the values of the nodes"""
        andNodes = model.andNodeList[
            currentNode
        ]  # find the list of shadow and nodes we must compute before computing value of current nodes
        andNodeInvertList = model.andNodeInvertList[
            currentNode
        ]  # find list of lists of whether input nodes need to be inverted (corresponds to inputOrder)
        writenode = (
            "" + model.nodeList[currentNode] + " , "
        )  # set up the initial string to use to write node
        inEdges = self.__findInEdges(model, currentNode)
        nodeIndividual = self.__simplifyRule(nodeIndividual, inEdges)
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
        """iterate over nodes to generate a BooleanNet representation for the entire model"""
        addString = ""
        for i in range(0, len(model.nodePositions)):
            addString = addString + model.writeNode_BoolNet(
                i,
                individual[model.individualParse[i] : model.individualParse[i + 1]],
                model,
            )
            addString = addString + "\n"
        return addString[:-1]

    def __eaMuPlusLambdaAdaptive(self, scSyncBoolC, graph, verbose=True):
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

        population = self.toolbox.population(n=self.params.popSize)

        logbook.header = ["gen", "nevals"] + (self.stats.fields if self.stats else [])
        lastcheck = []
        modellist = []
        fitnesslist = []
        popList = []

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in population if not ind.fitness.valid]
        # print("Invalid individuals: " + str(invalid_ind))
        updateBooler = ctypes.cdll.LoadLibrary("./simulator.so")
        scSyncBoolC = updateBooler.scSyncBool

        fitnesses = [
            indy[0]._ruleMaker__evaluateByNode(indy[1], KOlist, KIlist, scSyncBoolC)
            for indy in invalid_ind
        ]

        print("Fitnesses: " + str(fitnesses))
        print(len(fitnesses))

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
            return population, logbook

        # Begin the generational process
        for gen in range(1, ngen + 1):
            offspring = self.__varOrAdaptive(
                population, toolbox, lambda_, cxpb, mutpb, (1.0 * gen / ngen), mutModel
            )
            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            fitnesses = [
                indy[0]._ruleMaker__evaluateByNode(indy[1], KOlist, KIlist, scSyncBoolC)
                for indy in invalid_ind
            ]
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
                return population, logbook

        return population, logbook

    def __evaluateByNode(
        self,
        individual,
        KOlist,
        KIlist,
        cFunction,
        localSearch=False,
        importanceScores=False,
    ):
        """Includes Network Propagation"""
        model = self
        cellArray = []
        knockins = np.zeros(len(model.nodeList), dtype=np.intc, order="C")
        knockouts = np.zeros(len(model.nodeList), dtype=np.intc, order="C")

        for knocker in KOlist:
            knockouts[knocker] = 1
        for knocker in KOlist:
            knockins[knocker] = 1

        # put objects in correct format for passing to C
        nodeIndividual = np.array(individual, dtype=np.intc, order="C")
        indLen = len(nodeIndividual)
        andNodes = np.array(model.andNodes, dtype=np.intc, order="C")
        nodeNum = len(model.nodeList)
        andNodeInvert = np.array(model.andNodeInvert, dtype=object, order="C")
        individualParse = np.array(model.individualParse, dtype=np.intc, order="C")
        andLenList = np.array(model.andLenList, dtype=np.intc, order="C")
        nodePositions1 = model.nodePositions
        nodePositionsC = np.array(nodePositions1, dtype=np.intc, order="C")
        simSteps = self.params.simSteps
        lenSamples1 = len(model.sampleList)
        binMatC3 = np.array(
            copy.deepcopy(self.binMat.toarray(order="C")), order="C", dtype=np.intc
        )
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
            shape=(100, self.maxNodes), fill_value=0, dtype=np.intc, order="C"
        )  # simData[STEP][NODE]
        valsubmit = ctypes.c_void_p(vals.ctypes.data)
        lenSamples = ctypes.c_void_p(lenSamples1)

        localSearchC = ctypes.c_void_p(int(localSearch))
        importanceScoresC = ctypes.c_void_p(int(importanceScores))

        # errors = np.array(np.full(10000, fill_value=0, dtype=np.intc, order='C'))
        # errorsSubmit=ctypes.c_void_p(errors.ctypes.data)

        if localSearch:
            # look at errors node wise
            errors = np.array(
                np.full(self.maxNodes, fill_value=0, dtype=np.intc, order="C")
            )
            errorsSubmit = ctypes.c_void_p(errors.ctypes.data)
            cFunction(
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
                errorsSubmit,
                localSearchC,
                importanceScoresC,
            )  # in this case scSyncBoolC
            errors = errors.tolist()
            errors = errors[:nodeNum]
            return errors
        else:
            if importanceScores:
                importanceScores = np.array(
                    np.full(1, fill_value=0.0, dtype=np.float64, order="C")
                )
                importanceScoresC = ctypes.c_void_p(importanceScores.ctypes.data)
                cFunction(
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
                    importanceScoresC,
                )  # in this case importanceScore
                return importanceScores.tolist()
            else:
                # look at errors by sample
                errors = np.array(
                    np.full(self.maxSamples, fill_value=0, dtype=np.intc, order="C")
                )
                errorsSubmit = ctypes.c_void_p(errors.ctypes.data)
                cFunction(
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
                    errorsSubmit,
                    localSearchC,
                    importanceScoresC,
                )  # in this case scSyncBoolC
                errors = errors.tolist()
                return [sum(errors)]

    def __processERS(self, equivsName):
        """Create an individual from the ERS generated by the local search, for importance score calculation"""
        ersFile = open(str(equivsName), "rb")
        ers = pickle.load(ersFile)
        ersFile.close()

        # randomly sample the ers to make an individual
        individual = []
        for i in range(len(ers)):
            individual.extend(ers[i][randint(0, len(ers[i]) - 1)])
        return individual

    def __checkNodePossibilities(self, node, indy, KOlist, KIlist, scSyncBoolC):
        model = self
        tol = 0.0  # .01*len(newSSS) # set tolerance for equivalence
        end = model._ruleMaker__findEnd(node)  # find end of model for this node
        start = model.individualParse[node]  # find start of model for this node
        truth = list(indy[start:end])
        equivs = [truth]
        if (end - start) == 0:
            return truth, equivs, equivs, 0.0
        indOptions = []
        indErrors = []
        # iterate over possibilities for this node
        for i in range(1, 2 ** (end - start)):
            tempultimate = list(indy)
            tempInd = model._ruleMaker__bitList(i, len(truth))
            tempultimate[start:end] = tempInd  # set rule to one being checked
            currentsumtemp = self._ruleMaker__evaluateByNode(
                tempultimate, KOlist, KIlist, scSyncBoolC, localSearch=True
            )
            currentsum = currentsumtemp[
                node
            ]  # save the error found #subset complete error to
            indOptions.append(tempInd)
            indErrors.append(currentsum)
            gc.collect()
        minny = min(indErrors)
        equivs = []
        for i in range(len(indOptions)):
            if indErrors[i] <= minny + tol:  # changed from < to <=
                equivs.append(indOptions[i])
        truth = equivs[0]
        return (truth, equivs, minny, indErrors)

    def __calcImportance(self, equivs, model, importanceScore, graphName):
        # Create holder for importance scores
        importanceScoresDict = {}
        importanceScoreStdev = {}
        strat2_IS = {}
        strat3_IS = {}
        strat4_IS = {}
        tempList = list(range(0, len(self.nodeList)))
        # shuffle(tempList)
        # print(tempList)
        # print(len(tempList))
        for node in range(0, len(self.nodeList)):
            importanceScoresDict[self.nodeList[node]] = []
            importanceScoreStdev[self.nodeList[node]] = 0.0
        # Try 3 randomly sampled rule sets
        i = 0
        while i < 3:
            individual = self._ruleMaker__processERS(graphName + "_equivs1.pickle")
            for node in tempList:
                print(
                    "Node: "
                    + str(self.nodeList[node])
                    + ", Node Position: "
                    + str(node)
                )
                temp = self._ruleMaker__evaluateByNode(
                    individual,
                    [node],
                    [node],
                    importanceScore,
                    localSearch=False,
                    importanceScores=True,
                )
                print("Trial: " + str(i) + " Unprocessed IS: " + str(temp))
                importanceScoresDict[self.nodeList[node]].append(temp[0])
            i = i + 1

        print(importanceScoresDict)

        # Find maximum node importance score
        maxScore = max(importanceScoresDict.values())
        print("Max IS: " + str(maxScore))
        minScore = min(importanceScoresDict.values())
        print("Min IS: " + str(maxScore))

        # Rescaling to [0,1] using featureReScale
        for node in range(0, len(self.nodeList)):
            importanceScoresDict[self.nodeList[node]] = (
                importanceScoresDict[self.nodeList[node]][0] - minScore[0]
            ) / (maxScore[0] - minScore[0])

        print(importanceScoreStdev)

        ersFile = open(str(graphName + "_equivs1.pickle"), "rb")
        ers = pickle.load(ersFile)
        obsERS = {}
        maxERS = {}
        inDegreeNet = nx.read_graphml(graphName)
        # Normalize by number of rule sets that were tried
        for node in range(0, len(self.nodeList)):
            obsERS[self.nodeList[node]] = len(ers[node])
            inDegree = inDegreeNet.in_degree(self.nodeList[node])
            if inDegree == 0:
                maxERS[self.nodeList[node]] = 1
            else:
                inDegree = min(inDegree, 3)
                maxERS[self.nodeList[node]] = (
                    2 ** (len(ers[node][0])) - 1
                )  # 2**(inDegree+1) - 1 #
            # print(node)
            # print(obsERS[self.nodeList[node]])
            # print(maxERS[self.nodeList[node]])
            # Strategy 3: scale IS by (maxERS - obsERS + 1)/max ERS
            importanceScoresDict[self.nodeList[node]] = np.mean(
                importanceScoresDict[self.nodeList[node]]
            )
            importanceScoresDict[self.nodeList[node]] = importanceScoresDict[
                self.nodeList[node]
            ] * (
                (maxERS[self.nodeList[node]] - obsERS[self.nodeList[node]] + 1)
                / maxERS[self.nodeList[node]]
            )
            # importanceScoresDict[self.nodeList[node]] = np.mean(importanceScoresDict[self.nodeList[node]])
            # Strategy 1: divide by the standard deviation of the scores across rule set trials. This should be BEFORE rescaling to [0,1]
            # importanceScoreStdev[scObject.nodeList[node]] = np.std(importanceScoresDict[scObject.nodeList[node]])
            # newIS[scObject.nodeList[node]] = importanceScoresDict[scObject.nodeList[node]]/float(importanceScoreStdev[scObject.nodeList[node]])
            # Strategy 2: scale IS by log2((obsERS + 1)/max ERS)
            # strat2_IS[self.nodeList[node]] = importanceScoresDict[self.nodeList[node]] * (np.log2((obsERS[self.nodeList[node]] + 1)/maxERS[self.nodeList[node]]))
            # Strategy 3: scale IS by (maxERS - obsERS + 1)/max ERS
            # strat3_IS[self.nodeList[node]] = importanceScoresDict[self.nodeList[node]] * ((maxERS[self.nodeList[node]] - obsERS[self.nodeList[node]] + 1)/maxERS[self.nodeList[node]])
        # Strategy 4: abs(strat2)/max(strat2)
        # for node in range(0, len(self.nodeList)):
        #    strat4_IS[self.nodeList[node]] = abs(strat2_IS[self.nodeList[node]])/max([abs(val) for val in strat2_IS.values()])
        # Print out file of importance scores
        IS_df = pd.DataFrame(
            importanceScoresDict.items(), columns=["Node", "Importance Score"]
        )
        # IS_df["Std"] = IS_df["Node"].map(importanceScoreStdev)
        # IS_df["newScore"] = IS_df["Node"].map(newIS)
        # IS_df["Strat2_IS"] = IS_df["Node"].map#(strat2_IS)
        # IS_df["Strat3_IS"] = IS_df["Node"].map#(strat3_IS)
        # IS_df["Strat4_IS"] = IS_df["Node"].map(strat4_IS)
        IS_df["ObsERS"] = IS_df["Node"].map(obsERS)
        IS_df["MaxERS"] = IS_df["Node"].map(maxERS)

        IS_df.to_csv(
            str(graphName + "_importanceScores.csv"),
            sep=",",
            encoding="utf-8",
            index=False,
        )

        # Make graphml with importance scores as attributes
        net = self.ruleGraph
        nx.set_node_attributes(net, values=importanceScoresDict, name="importanceScore")
        # nx.set_node_attributes(net, values=strat2_IS, name='strat2_IS')
        # nx.set_node_attributes(net, values=strat3_IS, name='strat3_IS')
        # nx.set_node_attributes(net, values=strat4_IS, name='strat4_IS')
        nx.set_node_attributes(net, values=maxERS, name="maxERS")
        nx.set_node_attributes(net, values=obsERS, name="obsERS")

        # add abundance as attribute to graph
        binMat2 = self.binMat.A
        abundance = {}
        abundance_sd = {}
        numZeros = {}
        numOnes = {}

        for node in list(importanceScoresDict.keys()):
            node_index = self.geneList.index(node)
            expression = binMat2[node_index, :].tolist()
            abundance[node] = np.mean(expression)
            abundance_sd[node] = np.std(expression)
            expression = np.array(expression)
            numZeros[node] = (expression == 0).sum()
            numOnes[node] = (expression == 1).sum()

        nx.set_node_attributes(net, values=abundance, name="abundanceMean")
        nx.set_node_attributes(net, values=abundance_sd, name="abundanceStdev")
        nx.set_node_attributes(net, values=numZeros, name="abundanceZeros")
        nx.set_node_attributes(net, values=numOnes, name="abundanceOnes")

        nx.write_graphml_lxml(net, graphName[:-26] + "_IS.graphml")

        return importanceScoresDict

    class Params:
        def __init__(self):
            pass

        def __simParams(
            self,
            mutModel=0.25,
            cells=1,
            samples=1,
            generations=5,
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
