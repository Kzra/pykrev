from .calculate_mass import calculate_mass
from .element_counts import element_counts
import numpy as np
def kendrick_mass_defect(formula_list,mz_list,base = ['CH2']):
    """ 
	Docstring for function pyKrev.kendrick_mass_defect
	====================
	This function takes a list of molecular formula and their corresponding mz_values and calculates kendrick mass defect scores.  
    
	Use
	----
	approximate mass(Y,X)
    
	Returns a tuple contains two numpy arrays of len(Y). The first contains the kendrickMass and the second contains the kendrickMassDefect.
    
	Parameters
	----------
	Y: A list of molecular formula. 
    X: A list of mz_values of len(Y). 
    
    Base: Atom group used to define the Kendrick mass.
    Info
	----------
    Calculation taken from Hughey et al (2001) 
    "Kendrick Mass Defect Spectrum: A Compact Visual Analysis for Ultrahigh-Resolution Broadband Mass Spectra"
    
    """
    
    #Formula and mz lists must be the same length
    assert(len(formula_list)==len(mz_list))
    if type(mz_list) == list:
        mz_list = np.array(mz_list)
    if type(base) != list:
        raise AssertionError("Base must be provided as a string within a list. E.g. ['CH2']")

    kendrickMass = np.empty(len(formula_list))
    kendrickMassDefect = np.empty(len(formula_list))

    nominalBase = calculate_mass(base, method = 'nominal')
    exactBase = calculate_mass(base, method = 'monoisotopic')
    kendrickMass = np.array(mz_list) * (nominalBase/exactBase)
    kendrickMassDefect  = np.around(kendrickMass) - kendrickMass
    
    return kendrickMass, kendrickMassDefect