import numpy as np
from .element_counts import element_counts
def element_ratios(formula_list, ratios = ['OC','HC'], zero_ratio = True):
        
    """ 
	Docstring for function pyKrev.element_ratios
	====================
	This function takes a list of molecular formula strings and returns atomic ratios for each ratio given in ratios.
    
	Use
	----
	element_ratios(Y)
    
	Returns a list of len(Y) in which each item is a dictionary containing the atomic ratios listed in ratios, 
	calculated as the first element / second element. If the second element is zero, a nan or zerio value is returned
    depending on zero_ratio.
    
	Parameters
	----------
	Y: A list of molecular formula strings.
    
	ratios: A list of atomic ratios to compute, calculated as first element divided by second element. 
	Only C,H,N,O,P & S can be used. 
    
    zero_ratio: behaviour for if a zero ratio is accounted (i.e. the denominator is zero), if true a zero value is returned. 
    
    """
    count_list = element_counts(formula_list)
    ratio_list = []
    for count in count_list:
        ratio_counts = dict()
        for atom in ratios:
            if len(atom) > 2:
                raise Exception('Only C,H,N,O,P & S can be used. ')
            if count[atom[1]] == 0:
                if zero_ratio == True:
                    ratio_counts[atom] = 0
                else:
                    ratio_counts[atom] = np.nan
            else:
                ratio_counts[atom] = count[atom[0]]/count[atom[1]]
        ratio_list.append(ratio_counts)
        
    return ratio_list