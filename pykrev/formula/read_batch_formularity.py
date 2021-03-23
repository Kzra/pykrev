import pandas as pd
import numpy as np 
from ..diversity import ordination_matrix
def read_batch_formularity(report_name, ordination =True, impute_value = 'nan'):
    """ 
	Docstring for function PyKrev.read_batch_formularity
	====================
	This function reads in data from the report csv file produced by formularity software when performing batch assignment with alignment.
	The function filters out any formula that do not have C and H atoms.
    
	Use
	----
	read_batch_formularity(report_name)
    
	Returns a dataframe with the formula strings, masses, compound classes and mass errors of the overall alignment and their corresponding peak 
    intensities in each spectra, and an ordination matrix of type pd.Dataframe with formula strings as column headers and spectra as rows (see 
    pk.ordination_matrix) (optional) 
    
	Parameters
	----------
	report_name: name of csv file that the formularity report to be read is saved as. 
	impute_value: value to impute for 0 peak intensities in the output ordination matrix (default 0), see pk.Ordination_matrix)
    ordination: boolean, compute an ordination matrix using pk.ordination_matrix?
    """
    #read in csv to pandas
    report = pd.read_csv(report_name)
    #sum C and C13 to get total C number
    C = report['C'] + report['C13']
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
    
    ordinationMat = []
    if ordination == True:
        ordinationMat = ordination_matrix(molecular_formulas = spectraFormula,peak_intensities = spectraIntensities, group_names = spectraNames, impute_value = impute_value)
        #This will happen if there are identical formula in the molecular formula list
        if ordinationMat.shape[1] != formulaData.shape[0]:
            print('Warning, duplicate formula assignments detected. Ordination matrix will report values for first formula it encounters.')
    return formulaData, ordinationMat
