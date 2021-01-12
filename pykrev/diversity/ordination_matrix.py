import pandas as pd
import numpy as np
def ordination_matrix(molecular_formulas = [],peak_intensities=[],group_names = [], impute_value = 0):
    """ 
	Docstring for function pyKrev.ordination_matrix
	====================
	This function computes a sample data matrix for several lists of molecular formula nad corresponding peak intensities. 
    This matrix can be used for further ordination analysis (e.g. PCA, PcoA...)
    
	Use
	----
	ordination_matrix(molecular_formulas = [],peak_intensities =[],group_names = [])
    
	Returns a pandas dataframe, in which the column headers are a set of all formula across molecular_formula and the rows correspond to a specific sample.
    The [row,col] value of the dataframe is therefore the peak intensity of the formula in that sample. Impute value (default 0) if the formula was not present. 
    
	Parameters
	----------
	molecular_formulas: a list containing lists of molecular formulas
    peak_intensities: a list containing numpy arrays corresponding to the formula lists in molecular_formulas 
    group_names: a list of strings of len(molecular_formulas) corresponding to the row (i.e. sample) names
    impute_value: the value to impute when a formula isn't present in a group. An integer or float or 'nan' (default 0):

    
    """
    if  impute_value == 'nan':
        impute_value = np.nan
    if not group_names:
        group_names = [i for i in range(0,len(molecular_formulas))]
    all_formula = set()
    for flist in molecular_formulas:
        all_formula.update(flist) #update the set with the formula in molecular_formula
    all_formula = list(all_formula) #convert back to list for ease with indexing 
    ordination_mat = pd.DataFrame(columns = all_formula,index = group_names)
    col_index = 0
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