import pandas as pd
import numpy as np 
from .msTuple import msTuple
from .msTupleDict import msTupleDict
from ..diversity.ordination_matrix import ordination_matrix
def read_batch_formularity(report_name):
    """ 
	Docstring for function PyKrev.read_batch_formularity
	====================
	This function reads in data from the report csv file produced by formularity software when performing batch assignment with alignment.
	The function filters out any formula that do not have C and H atoms.
    
	Use
	----
	read_batch_formularity(report_name)
    
	Returns an msTupleDict contaning sample names as keys and corresponding msTuples as values 

    
	Parameters
	----------
	report_name: name of csv file that the formularity report to be read is saved as. 
    
    Info
    -----------
    PyKrev will filter out formula with 13C assignments
    """
    #Main
    report = pd.read_csv(report_name)
    notIsotopologue = report['C13'] == 0 # boolean array of only non isotopologues
    report = report[notIsotopologue]
    report.reset_index(drop = True, inplace = True) #reset the index to account for the removed isotopologues]
    C = report['C']
    H = report['H']
    O = report['O']
    N = report['N'] 
    S = report['S']
    P = report['P']
    #dimensions of report (number of formula, samples)
    x,y = report.shape
    #list to store all formula in
    molecular_formula = []
    #turn atom numbers into formula strings 
    for n in range(0,len(C)):
        if C[n] == 0 and H[n] == 0: #If there isn't a count for C and H don't include the formula
            temp_formula = np.nan
        else:
            temp_formula = 'C'+str(C[n])+'H'+str(H[n]) #add C and H then ... 
            if N[n] > 0:
                temp_formula = temp_formula + 'N' + str(N[n])
            if O[n] > 0:
                temp_formula = temp_formula + 'O' + str(O[n])
            if P[n] > 0:
                temp_formula = temp_formula + 'P' + str(P[n])
            if S[n] > 0:
                temp_formula = temp_formula + 'S' + str(S[n])
        molecular_formula.append(temp_formula)
    #convert to numpy array so we can use boolean indexing
    molecular_formula = np.array(molecular_formula)
    #remove nan values
    filter_formula = molecular_formula[molecular_formula != 'nan']
    #add to the original dataframe to preserve indexing
    report['all_formula'] = molecular_formula
    #make new dataframe
    formulaData = pd.DataFrame(index=None)
    #copy across formula, mass, cclass and mass error data with boolean indexing
    formulaData['formula'] = report['all_formula'][molecular_formula != 'nan']
    formulaData['mass'] = report['Mass'][molecular_formula != 'nan']
    formulaData['cclass'] = report['Class'][molecular_formula != 'nan']
    formulaData['massError'] = report['Error_ppm'][molecular_formula != 'nan']
    spectraNames = []
    spectraFormula = []
    spectraIntensities = []
    #column 15 is the first sample ID column
    for i in range(14,y):
        #save the peak intensities after boolean indexing has been applied
        intensity = report.iloc[:,i][molecular_formula != 'nan']
        formulaData[report.columns[i]] = intensity
        spectraNames.append(report.columns[i])
        spectraIntensities.append(np.array(intensity[intensity > 0]))
        spectraFormula.append(list(filter_formula[intensity > 0]))
    #Now process the formulaData into an msTuple dictionary 
    batchDict = msTupleDict()
    x,y = formulaData.shape
    sampleNames = formulaData.columns[4::]
    for name in sampleNames:
        batchDict[name] = msTuple(list(formulaData['formula'][formulaData[name] > 0]), np.array(formulaData[name][formulaData[name] > 0]), np.array(formulaData['mass'][formulaData[name] > 0]))
    return batchDict
