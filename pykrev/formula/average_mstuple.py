import numpy as np
from .calculate_mass import calculate_mass
from .msTuple import msTuple
def average_mstuple(Y, intensityMethod = 'mean', mzMethod = 'mean', minOccurrence = 1, zeroValues = True, stdDev = False):
    """ 
    Docstring for function pykrev.average_mstuple
    ==========
    Averages all the spectra in an msTupleDict and returns the average msTuple
    
    Use 
    ----------
    average_mstuple(Y) 
    
    Returns an msTuple with the spectra in the msTuple dictionary averaged according to their intensity and mz. 
    
    Parameters 
    ----------
    Y: An msTupleDict
    
    Give a description of all possible input parameters and their type, starting with positional arguments. 

    intensityMethod: string, one of: 
        'mean' - use the mean intensity for a given formula 
        'max' - use the max intensity for a given formula 
        'sum' - use the sum intensity for a given forumla
    
    mzMethod: string, one of:
        'mean' - use the mean mz for each formula 
        'median' - use the median mz for each formula
        'monoisotopic' - use the monoisotopic mz for each formula 

    minOccurance: int, the minimum number of samples a formula must appear in to be included in the intensity averaging

    zeroValues: bool, include zero (i.e. missing formula) values in average intensity calculation

    stdDev: boolean, if true return also a tuple containing 
                    (i)  array of the sample tandard deviations of each formula intensity in the output msTuple, 
                    (ii) array of the standard deviatons of each mz in the output msTuple

    """
    #Tests
    Y.validate()
    assert intensityMethod in ['mean','max','sum'], "You must provide a valid method"
    assert mzMethod in ['mean','median', 'monoisotopic'], "You must provide a valid method"
    assert minOccurrence > 0, "minOccurrence must be at least 1"
    #Setup
    OrdinationMat = Y.to_OrdinationMatrix(impute_value = 0)
    row,col = OrdinationMat.shape
    formulaNames = OrdinationMat.columns
    outputFormula = []
    outputIntensity = []
    outputMZ  = []
    stdDevIntensity = []
    stdDevMZ = []
    #Main
    for i in range(col):
        ##Test whether Occurrence matches min Occurrence
        nonZeroIdx = OrdinationMat.iloc[:,i].to_numpy().nonzero()[0]
        occurrence = len(nonZeroIdx)
        if occurrence >= minOccurrence:
            ## Remove or keep zero values
            if zeroValues == True:
                testArray = OrdinationMat.iloc[:,i].to_numpy()
            elif zeroValues == False:
                testArray = OrdinationMat.iloc[:,i].to_numpy()[nonZeroIdx]
            ## Perform the averaging
            if intensityMethod == 'mean':
                outputIntensity.append(testArray.mean())
            elif intensityMethod == 'max':
                outputIntensity.append(testArray.max())                
            elif intensityMethod == 'sum':
                outputIntensity.append(testArray.sum())
            formula = formulaNames[i]
            outputFormula.append(formula)
            ## Determine the mz across spectra where this formula was found 
            if mzMethod == 'monoisotopic':
                outputMZ.append(calculate_mass(formula))
            mzArray = np.array([])
            for v in Y.values():
                try:
                    idx = v.formula.index(formula)
                except ValueError:
                    continue
                mzArray = np.append(mzArray,v.mz[idx])
            if mzMethod == 'mean':
                outputMZ.append(mzArray.mean())
            elif mzMethod == 'median':
                outputMZ.append(mzArray.median())
            ## Determine the standard deviations
            if stdDev == True:
                stdDevIntensity.append(testArray.std(ddof = 1))
                stdDevMZ.append(mzArray.std(ddof = 1))
    assert len(outputIntensity) == len(outputMZ) == len(outputFormula)
    outputMZ = np.array(outputMZ)
    outputIntensity = np.array(outputIntensity)
    if stdDev == True:
        assert len(outputIntensity) == len(stdDevIntensity)
        return msTuple(outputFormula, outputIntensity, outputMZ), (stdDevIntensity, stdDevMZ)
    else:
        return msTuple(outputFormula, outputIntensity, outputMZ)




            



