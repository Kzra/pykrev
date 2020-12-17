import numpy as np
def filter_spectral_interference(mass_list,formula_list,peak_intensities,tol = 2,verbose = True):
    """ 
	Docstring for function pykrev.filter_spectral_interference
	====================
	This function takes a list or numpy array of masses and a corresponding list of formula strings 
    and a corresponding list or numpy array of peak intensities
    and returns filtered copies of all arrays with suspected 
    spectral interferences due to high molecular weight doubly charged ions, removed.
    
	Use
	----
	filter_spectral_interference(Y,X,Z)
    
	Returns a tuple containing an array of len(Y[filter]) in which each item is the filtered mass, and a list of len(Y[filter]) in which each item 
    is the filtered formula and an array of len(Y[filter]) in which each item is a filtered peak_intensity. 
    
	Parameters
	----------
	Y: A list or numpy array of masses.
    X: A list of molecular formula strings.
    Z: A list or numpy array of peak intensities.
    tol: Integer of float,the tolerance (in ppm) used to identify monoisotopic 
    verbose: Booloean, print text to the console.
    
    Note: for more information on this algorithm please refer to Patriarca, Claudia, and Jeffrey A. Hawkes. 
    "High Molecular Weight Spectral Interferences in Mass Spectra of Dissolved Organic Matter." 
    Journal of the American Society for Mass Spectrometry (2020).
    """    
    assert len(mass_list) == len(formula_list) == len(peak_intensities), "All input arrays must be of equal size"
    mass_list = np.array(mass_list)
    formula_list = np.array(formula_list)
    peak_intensities = np.array(peak_intensities)
    mass_defect = mass_list - np.floor(mass_list)
    c13d2 = (13.00335-12)/2
    flags = np.zeros(len(mass_list))
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
        print(f"{len(mass_list)-len(filter_mass)} masses removed.")
    
    return filter_mass, filter_formula, filter_peak_intensities