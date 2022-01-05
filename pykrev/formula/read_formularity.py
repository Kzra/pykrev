import pandas as pd
import numpy as np 
from .msTuple import msTuple
def read_formularity(report_name,pi_col = []):
    """ 
	Docstring for function PyKrev.read_formularity
	====================
	This function reads the report csv file produced by formularity software and returns an msTuple.
    
	Use
	----
	read_formularity(report_name)
    
	Returns an msTuple. 

	Parameters
	----------
	report_name: name of csv file that the formularity report to be read is saved as. 
	pi_col: name of the column in that file that peak intensities are found in. If not given the last column is used. 

    Note: PyKrev will filter out formula with 13C assignments
    """
    report = pd.read_csv(report_name)
    notIsotopologue = report['C13'] == 0 # boolean array of only non isotopologues
    report = report[notIsotopologue]
    report.reset_index(drop = True, inplace = True) #reset the index to account for the removed isotopologues
    C = report['C']
    H = report['H']
    O = report['O']
    N = report['N'] 
    S = report['S']
    P = report['P']
    if not pi_col: 
        I = report.iloc[:,-1] #take the final column of the report file to contain peak intensities. Not sure how stable this is.
    else:
        I = report[pi_col] #user supplied a column name for the peak intensities.
    MZ = report['Mass']
    G = report['Class']
    ME = report['Error_ppm']
    molecular_formula = []
    mass_charge = np.array([])
    mass_error = np.array([])
    peak_intensities = np.array([])
    compound_class = []
    for n in range(0,len(C)):
        if C[n] == 0 and H[n] == 0: #If there isn't a count for C and H don't include the formula
            pass
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
            peak_intensities = np.append(peak_intensities,I[n])
            mass_charge =  np.append(mass_charge,MZ[n])
            mass_error = np.append(mass_error,ME[n])
            compound_class.append(G[n])
    return msTuple(molecular_formula, peak_intensities, mass_charge) #return an msTuple