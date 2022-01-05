from ..formula.calculate_mass import calculate_mass
import networkx as nx
import numpy as np
def reaction_network(msTuple, filePath = '', fileFormat = 'none', reactionDict = {
                                            'decarboxylation': -calculate_mass(['CO2']),
                                            'methylation': calculate_mass(['CH2']),
                                            'demethylation': -calculate_mass(['CH2']),
                                            'hydrogenation': calculate_mass(['H2']),
                                            'dehydrogenation': -calculate_mass(['H2']),
                                            'hydration': calculate_mass(['H2O']),
                                            'dehydration': -calculate_mass(['H2O']),
                                            'oxidation': calculate_mass(['O']),
                                            'reduction': -calculate_mass(['O'])
                                            }, nodeAnnotations = {}, roundVal = 8):
    """ 
	Docstring for function PyKrev.reaction_network
	====================
	This function takes an msTuple and writes a directed graph format network representation to the location set by filePath using the networkx library. 
    The reaction network has molecular formula as nodes and the reactions in reactionDict as edges. 
    The nodes can be annotated by the user with custom values.
    Edges are annotated with the associated reaction name.
    
	Use
	----
	reactionGraph(Y)
    
	Returns a tuple containing (i) a networkx representation of the graph and (ii) a dictionary of reaction counts 
    
	Parameters
	----------
	Y: msTuple
    reactionDict: dictionary, containing reaction names as keys and their associated change in monoisotopic formula mass as values.
    nodeAnnotations: dictionary, containing the node annotation names as keys, and their associated values as numpy arrays. e.g. {'Peak Intensity': intensityArray)
        where len(intensityArray) == len(formulaList)
    filePath: string, directory location to write the graph file to.
    fileFormat: string, the format of the graph file to write, either 'grapml' or 'gexf' or 'none'. If 'none' no graph file is written.
    roundVal: int, number of decimal places to round to in the mass defect calculation, default is 8
    """ 
    #Tests
    assert fileFormat in ['graphml', 'gexf', 'none'], "format must be graphml or gexf"
    assert type(filePath) == str, "filePath must be provided as a string"
    formulaList = msTuple[0]
    for value in nodeAnnotations.values():
        assert len(value) == len(formulaList), "ensure the arrays in nodeAnnotations are the same length as formula list"
    #Setup
    graphNodes = []
    reactionList = list(reactionDict.keys())
    reactionCounts = reactionDict.fromkeys(reactionDict, 0)
    formulaMass = calculate_mass(formulaList)
    #Main
    for i in range(0,len(formulaList)):
        # first process the node annotations into a usable format 
        nodeDict = {}
        for key in nodeAnnotations.keys():
            nodeDict[key] = nodeAnnotations[key][i]
        graphNodes.append((formulaList[i],nodeDict)) #set up the graph nodes with nodeAnnotations 
    G = nx.DiGraph() # Use a directed graph (edges have directional information)
    G.add_nodes_from(graphNodes) #Each formula in our sample is a node
    for j in range(0,len(formulaList)):
        for reactionType in reactionList:
            matchIndex = np.nonzero(np.round(formulaMass, roundVal) == np.round((formulaMass[j] + reactionDict[reactionType]), roundVal))
            if matchIndex[0].size > 0:
                reactionCounts[reactionType] += 1
                G.add_edge(formulaList[j],formulaList[int(matchIndex[0])], reaction = reactionType) 
    print('Number of nodes in graph:')
    print(G.number_of_nodes())
    print('Number of edges in graph')
    print(G.number_of_edges())
    if fileFormat == 'graphml':
        nx.write_graphml(G, filePath)
    elif fileFormat == 'gexf':
        nx.write_gexf(G, filePath)
    return G, reactionCounts
