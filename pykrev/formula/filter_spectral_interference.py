import numpy as np
from .msTuple import msTuple
def filter_spectral_interference(msTupleObj, tol = 2 , verbose = True):
    """ 
	Docstring for function pykrev.filter_spectral_interference
	====================
	This function takes an msTuple and returns filtered copies of all arrays with suspected spectral interferences due to high molecular weight doubly charged ions, removed.
    
	Use
	----
	filter_spectral_interference(Y)
    
	Returns an msTuple containing filtered arrays.
    
	Parameters
	----------
    Y: msTuple
    tol: Integer of float,the tolerance (in ppm) used to identify monoisotopic 
    verbose: Booloean, print text to the console.

    Info
    -----------
    For more information on this algorithm refer to Patriarca, Claudia, and Jeffrey A. Hawkes. 
    "High Molecular Weight Spectral Interferences in Mass Spectra of Dissolved Organic Matter." 
    Journal of the American Society for Mass Spectrometry (2020).
    """    
    #Setup
    mass_list = msTupleObj[2]
    formula_list = np.array(msTupleObj[0])
    peak_intensities = msTupleObj[1]
    mass_defect = mass_list - np.floor(mass_list)
    c13d2 = (13.00335-12)/2
    flags = np.zeros(len(mass_list))
    #Main
    for i in range(0,len(mass_list)):
        if mass_defect[i] > 0.4 and mass_defect[i] < 0.8:
            flags[i] = 1
            mono_mass = mass_list[i] - c13d2
            mass_diff = abs(mass_list - mono_mass)
            mono_pos = mass_diff == min(mass_diff)
            if min(mass_diff)/mass_list[i] * 1e6 < tol:
                if peak_intensities[mono_pos]/peak_intensities[i] < 10:
                    flags[mono_pos] = 2
    spectralFilter = flags == 0 
    filter_mass = mass_list[spectralFilter]
    filter_formula = formula_list[spectralFilter]
    filter_peak_intensities = peak_intensities[spectralFilter]
    if verbose == True:
        print(f"{len(mass_list)-len(filter_mass)} interferences removed.")
    filter_formula = list(filter_formula) #reconvert filter_formula to list to sort out data type issue
    return msTuple(filter_formula, filter_peak_intensities, filter_mass)