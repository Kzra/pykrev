import numpy as np
def nominal_oxidation_state(count_list):
    
    """ 
	Docstring for function pyKrev.nominal_oxidation_state
	====================
	This function takes a list of element counts and gives the nominal oxidatate state of C.
    
	Use
	----
	nominal_oxidation_state(Y)
    
	Returns a numpy array of len(Y) in which each item is the nominal oxidation state of C.  
    
	Parameters
	----------
	Y: A list of dictionary items containing atomic counts as produced by element_counts. 
	Must include C,H,O,N,P,S. 
    
    
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
    assert isinstance(count_list,list),'supply a list of counts given by element_counts()'
    assert all(isinstance(i,dict) for i in count_list), 'supply a list of counts given by element_counts()'
    
    NOSCs = np.array([])
        
    for i in count_list:
        NOSCs = np.append(NOSCs,-((4*i['C'] + i['H'] - 3 * i['N'] - 2 * i['O'] + 5 * i['P'] - 2 * i['S'])/i['C']) + 4) 
    #still need to test this calculation is correct
    
    return NOSCs