import pandas as pd
import numpy as np
def ordination_matrix(msTupleDict, impute_value = 'nan'):
    """ 
	Docstring for function pyKrev.ordination_matrix
	====================
	This function computes a sample data matrix from an msTupleDict
    This matrix can be used for further ordination analysis (e.g. PCA, PCoA...)
    
	Use
	----
	ordination_matrix(Y)
    
	Returns a pandas dataframe in which the column headers are a set of all formula found in msTupleDict and the rows correspond to a specific sample.
    The [row,col] value of the dataframe is therefore the peak intensity of a formula. Impute value (default 0) if the formula was not present. 
    
	Parameters
	----------
	Y: an msTupleDict
    impute_value: the value to impute when a formula isn't present in a group. An integer or float or 'nan' (default 0):
    """
    #Setup
    if  impute_value == 'nan':
        impute_value = np.nan
    group_names = msTupleDict.keys()
    molecular_formulas = []
    peak_intensities = []
    for msTuple in msTupleDict.values():
        molecular_formulas.append(msTuple[0])
        peak_intensities.append(msTuple[1])
    all_formula = set()
    for flist in molecular_formulas:
        all_formula.update(flist) #update the set with the formula in molecular_formula
    all_formula = list(all_formula) #convert back to list for ease with indexing 
    ordination_mat = pd.DataFrame(columns = all_formula,index = group_names)
    col_index = 0
    #Main
    for f in all_formula: #cycle through all_formula 
        row_index = 0 
        for flist, plist in zip(molecular_formulas,peak_intensities):
            try:
                i = flist.index(f) #find the first instance of the formula in flist
                ordination_mat.iloc[row_index,col_index] = plist[i] #set the row value to the peak intensity for that sample 
                row_index += 1  #move down a row i.e. sample
            except ValueError: #if the formula can't be found 
                ordination_mat.iloc[row_index,col_index] = impute_value #set the value to impute_value 
                row_index += 1 #move down a row
        col_index += 1 #move across to the next column i.e. formula
    return ordination_mat
