import numpy as np
from .element_counts import element_counts
def calculate_mass(msTuple, method = 'monoisotopic', protonated = False, ion_charge = 0):     
    """ 
	Docstring for function pyKrev.calculate_mass
	==========
	This function takes an msTuple and calculates the monoisotopic, nominal or average mass of each formula in the formula list.
    
	Use
	----------
	calculate_mass(Y)
    
	Returns a numpy array of len(Y[0]) in which each item , i , is the calculated mass of Y[0][i].
    
	Parameters
	----------
	Y: msTuple OR a list of molecular formula strings
    
    Method: String, the type of mass calculation to perform. One of:
        - 'monoisotopic', i.e. based on the exact mass of the most abundant isotope. 
        - 'nominal', i.e. the mass of the most abundant isotope rounded to the nearest integer.
        - 'average', i.e. mass taking into account weighted abundance of all natural isotopes of the element.

    protonated: Boolean, if True calculate mass for close shell ions([M + H]+  or [M - H]-) depending on ion_charge

    ion_charge: int, the ion charge of the formula, if not charged == 0 

	Info
	----------
    average masses taken from: http://physics.nist.gov/PhysRefData/Compositions/.
    monoisotopic masses taken from CoreMS python library.
    """
    #Setup 
    count_list = element_counts(msTuple)
    electron_mass = 0.0005485
    if method == 'monoisotopic':
        element_masses = {'C': 12.0, 'H': 1.007825032239, 'O':15.9949146195717,'N': 14.003074004432,'S': 31.972071174414, 'P': 30.973761998427, 'Cl': 34.96885268237, 'F': 18.9984031627392}
    elif method == 'nominal':
        element_masses = {'C': 12, 'H': 1, 'O':16,'N': 14,'S': 32, 'P': 31, 'Cl':35, 'F':19}
    elif method == 'average':
        element_masses = {'C': 12.010736, 'H': 1.007941, 'O':15.999405,'N': 14.006743,'S': 32.066085, 'P': 30.973762, 'Cl': 35.452938, 'F': 18.998403}
    else:
        raise Exception("Method not recognised.")
    mass_list = np.empty(len(count_list))
    #Main
    for i in range(0,len(count_list)):
        mass_list[i] = count_list[i]['C'] * element_masses['C'] + count_list[i]['H'] * element_masses['H'] + count_list[i]['N'] * element_masses['N'] + count_list[i]['O'] * element_masses['O'] + count_list[i]['P'] * element_masses['P'] + count_list[i]['S'] * element_masses['S'] + count_list[i]['Cl'] * element_masses['Cl'] + count_list[i]['F'] * element_masses['F']
        if ion_charge != 0: 
            mass_list[i] += (ion_charge * -1 * electron_mass)
            if protonated == True:
                mass_list[i] += (ion_charge * element_masses['H'])
            mass_list[i] = mass_list[i]/abs(ion_charge)
    return mass_list
