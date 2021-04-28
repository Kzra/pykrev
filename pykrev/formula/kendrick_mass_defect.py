from .calculate_mass import calculate_mass
from .element_counts import element_counts
import numpy as np
def kendrick_mass_defect(formula_list,mz_list,base = ['CH2'], rounding = 'even'):
    """ 
	Docstring for function pyKrev.kendrick_mass_defect
	====================
	This function takes a list of molecular formula and their corresponding mz_values and calculates kendrick mass defect scores.  
    
	Use
	----
	kendrick_mass_defect(Y,X)
    
	Returns a tuple contains two numpy arrays of len(Y). The first contains the kendrickMass and the second contains the kendrickMassDefect.
    
	Parameters
	----------
	Y: A list of molecular formula. 
    X: A list of mz_values of len(Y). 
    
    base: Atom group used to define the Kendrick mass.
    rounding: The method of rounding to use when calculating kendrick mass defect. One of 
        - even: for values exactly halfway between rounded decimal values, numpy rounds to the nearest even value. 
                Thus 1.5 and 2.5 round to 2.0, -0.5 and 0.5 round to 0.0, etc.
        - integer: for values exactly halfway between rounded decimal values, numpy rounds up to the nearest integer value. 
                Thus 1.5 becomes 2.0, 2.5 becomes 3.0 
    Info
	----------
    Calculation taken from Hughey et al (2001) 
    "Kendrick Mass Defect Spectrum: A Compact Visual Analysis for Ultrahigh-Resolution Broadband Mass Spectra"
    Note: Rounding calculations may lead to "edge effects" in complex datasets with many peaks. 
    We recommend experimenting with both even and integer rounding methods when making kmd plots.

    """
    
    #Formula and mz lists must be the same length
    assert(len(formula_list)==len(mz_list))
    if type(mz_list) == list:
        mz_list = np.array(mz_list)
    if type(base) != list:
        raise AssertionError("Base must be provided as a string within a list. E.g. ['CH2']")
    #Check for valid rounding method 
    assert(rounding in ['even','integer']), "provide a valid rounding method"

    kendrickMass = np.empty(len(formula_list))
    kendrickMassDefect = np.empty(len(formula_list))

    nominalBase = calculate_mass(base, method = 'nominal')
    exactBase = calculate_mass(base, method = 'monoisotopic')
    kendrickMass = np.array(mz_list) * (nominalBase/exactBase)
    if rounding == 'even':
        kendrickMassDefect  = np.around(kendrickMass) - kendrickMass
    elif rounding == 'integer':
        kendrickMassDefect  = np.rint(kendrickMass) - kendrickMass
    
    return kendrickMass, kendrickMassDefect