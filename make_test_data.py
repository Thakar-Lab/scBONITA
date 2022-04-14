# process scbonita test data
seed(15)
testData = {}
testData["originalGraph"] = nx.read_graphml("hsa00010.graphml")
testData["processedGraph"] = nx.read_graphml("hsa00010.graphml_processed.graphml")
testData["trainingData"] = pd.read_csv("subpop_14.bin", index_col=0)
testData["importanceScores"] = pd.read_csv(
    "hsa00010.graphml_processed.graphml_importanceScores.csv"
)
testData["scObject"] = pickle.load(open("subpop_14.binscTest.pickle", "rb"))
testData["metaData"] = pd.read_csv("allMetaData.txt", sep=",", index_col=0)
testData["contrast"] = pd.read_csv("contrasts.txt", sep="\t", header=None)
testData["embeddings"] = pd.read_csv("cell_embeddings.csv", index_col=0)
# sample dataset to create a dataset of 100 cells
originalColumnNames = sample(list(testData["trainingData"].columns), 100)
temp = testData["trainingData"].loc[:, originalColumnNames]
# rename cells
temp.columns = ["Cell" + str(i) for i in list(range(0, len(temp.columns)))]
temp.to_csv("trainingData.csv")
originalColumnNames
testData["metaData"].loc[list(originalColumnNames), "PlaqueScore"]
tempConditions = pd.DataFrame(
    testData["metaData"].loc[list(originalColumnNames), "PlaqueScore"]
)
tempConditions.index = [
    "Cell" + str(i) for i in list(range(0, len(temp.columns)))
]  # originalColumnNames
# revalue conditions.loc
tempConditions.columns = ["Conditions"]
tempConditions.Conditions = [
    "Treatment" if i == "AS+" else "Control" for i in tempConditions.Conditions
]
tempConditions.to_csv("conditions.txt")
tempEmbeddings = testData["embeddings"].loc[list(originalColumnNames)]
tempEmbeddings.index = ["Cell" + str(i) for i in list(range(0, len(temp.columns)))]
tempEmbeddings
tempEmbeddings.to_csv("testData_embeddings.csv")
