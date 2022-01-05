import numpy as np
from .element_counts import element_counts
def nominal_oxidation_state(msTuple):
    """ 
	Docstring for function pyKrev.nominal_oxidation_state
	====================
	This function takes an msTuple and returns the nominal oxidatate state of C.
    
	Use
	----
	nominal_oxidation_state(Y)
    
	Returns a numpy array of len(Y[0]) in which each item, i, is the nominal oxidation state of C corresponding to Y[0][i].  
    
	Parameters
	----------
	Y: msTuple
    
	Info
	----------
	This function derives the nominal oxidation state of C as described in LaRowe, Douglas E., and Philippe Van Cappellen. 
	"NOSC allows for computation of the average carbon oxidation state of an organic compound without writing out the oxidation half reaction."
	"Degradation of natural organic matter: a thermodynamic analysis." 
	Geochimica et Cosmochimica Acta 75.8 (2011): 2030-2042.       

	NOSC = -((-Z + 4a + b - 3c - 2d + 5e - 2f)/a) + 4 
	where Z = charge and 
	a = C
	b = H
	c = N 
	d = O
	e = P
 	f = S
        
    """      
	#Setup
    count_list = element_counts(msTuple)
    NOSCs = np.array([])
    #Main
    for i in count_list:
        NOSCs = np.append(NOSCs,-((4*i['C'] + i['H'] - 3 * i['N'] - 2 * i['O'] + 5 * i['P'] - 2 * i['S'])/i['C']) + 4) 
    return NOSCs