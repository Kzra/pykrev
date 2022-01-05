from ..formula.calculate_mass import calculate_mass
import numpy.linalg as la
import numpy as np
def page_rank(msTuple, reactionDict = {
                                            'decarboxylation': -calculate_mass(['CO2']),
                                            'methylation': calculate_mass(['CH2']),
                                            'demethylation': -calculate_mass(['CH2']),
                                            'hydrogenation': calculate_mass(['H2']),
                                            'dehydrogenation': -calculate_mass(['H2']),
                                            'hydration': calculate_mass(['H2O']),
                                            'dehydration': -calculate_mass(['H2O']),
                                            'oxidation': calculate_mass(['O']),
                                            'reduction': -calculate_mass(['O'])
                                            }, reactionWeights = {}, d = 0.9, tol = 0.01, roundVal = 8):
    """ 
	Docstring for function PyKrev.page_rank
	====================
	This function takes an msTuple and performs the pagerank algorithm on a reaction network derived from  the list of formula. 
    The reaction network can have different weights for each reaction type this should be provided in the dictionary reaction weights.
    For more information on the page rank algorithm and it's implementation here, see the docs: pykrev/docs/pagerank_and_networkvis/PageRankandNetworkVis.ipynb
    
	Use
	----
	page_rank(Y)
    
	Returns a numpy array of len(Y[0]) with the pagerank scores corresponding to the elements in Y[0].
    
	Parameters
	----------
	Y: msTuple
    reactionDict: dictionary, containing reaction names as keys and their associated change in monoisotopic formula mass as values.
    reactionWeights: dictionary, containing the relative weighting to give to each reactionType. If not provided each reactionWeight is given with equal value.
    d: float, damping factor in page rank algorithm
    tol: float, tolerance to run power iteration method to
    roundVal: int, number of digits to round to gor mass defect calculations
    """ 
    #Tests
    if len(reactionWeights) == 0:
        for key in reactionDict.keys():
            reactionWeights[key] = 1
    else: 
        assert reactionWeights.keys() == reactionDict.keys(), "reactionWeights and reactionKeys must have identical keys"
    #Setup
    formulaList = msTuple[0]
    N = len(formulaList)
    formulaMass = calculate_mass(formulaList) # Compute the exact monoisotopic mass of each formula in the dataset
    reactionList = list(reactionDict.keys())
    #Main
    ## Create the matrix L 
    L = np.zeros([N,N])
    for j in range(N):
        for reactionType in reactionList:
            # find possible matches, use np.round to account for rounding errors
            matchIndex = np.where(np.round(formulaMass,roundVal) == np.round((formulaMass[j] + reactionDict[reactionType]),roundVal))
            # if a match is found
            if len(matchIndex[0]) > 0:
                L[matchIndex,j] = reactionWeights[reactionType]
        ## normalise the probabilities so they sum to one
        if sum(L[:,j]) > 0:
            L[:,j] = L[:,j]/sum(L[:,j])
        ## if the column sums to zero set to 1/N
        else: 
            L[:,j] = 1/N
    # run the pageRank algorithm to convergence
    r = 100 * np.ones(N) / N # Sets up the probability vector
    M = d * L + (1-d)/N * np.ones([N, N]) 
    lastR = r
    r = M @ r
    i = 0
    while la.norm(lastR - r) > tol :
        lastR = r
        r = M @ r
        i += 1
    print(str(i) + " iterations to convergence.")
    return r