from .element_counts import element_counts
import numpy as np
def double_bond_equivalent(msTuple):
    """ 
	Docstring for function pyKrev.double_bond_equivalent
	==========
	This function takes an msTuple and returns the double bond equivalent of each formula in the formula list.
    
	Use
	----------
	double_bond_equivalent(Y)
    
	Returns a numpy array of len(Y[0]) in which each item , i , is the double bond equivalent of Y[0][i].  
    
	Parameters
	----------
	Y: msTuple OR a list of molecular formula strings
    
	Info
	----------
	Double bond equivalent (DBE; UN; degree of unsaturation; PBoR [Pi Bonds or Rings]): 
	The number of molecules of H2 that would have to be added to a molecule to convert all pi bonds to single bonds, 
	and all rings to acyclic structures. 
	The DBE number can be calculated from the formula using the following equation:
	DBE = UN = PBoR = C - (H/2) + (N/2) +1,
	where: C = number of carbon atoms, H = number of hydrogen and halogen atoms, and N = number of nitrogen atoms.
    """    
    #Setup
    count_list = element_counts(msTuple)
    DBE_array = np.array([])
    #Main
    for count in count_list:
        Halogens = ['H','Cl','Br','I','F','At','Ts']
        Hal = 0 
        for el in Halogens: 
            try:
                Hal += count[el]
            except KeyError:
                pass 
            
        DBE_counts = count['C'] - (Hal/2) + (count['N']/2) + 1 
        DBE_array = np.append(DBE_array,DBE_counts)        
    return DBE_array