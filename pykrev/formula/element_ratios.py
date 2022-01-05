import numpy as np
from .element_counts import element_counts
def element_ratios(msTuple, ratios = ['OC','HC'], zero_ratio = True):
        
    """ 
	Docstring for function pyKrev.element_ratios
	====================
	This function takes an msTuple and returns atomic ratios of each ratio given in ratios for each formula in the formula list.
    
	Use
	----
	element_ratios(Y)
    
	Returns a list of len(Y) in which each item is a dictionary containing the atomic ratios listed in ratios, 
	calculated as the first element / second element. If the second element is zero, a nan or zero value is returned
    depending on zero_ratio.
    
	Parameters
	----------
	Y: msTuple OR a list of molecular formula strings
    
	ratios: List, atomic ratios to compute, calculated as first element divided by second element. Only C,H,N,O,P & S can be used. 
    
    zero_ratio: Boolean, behaviour for if a zero ratio is encountered (i.e. the denominator is zero), if true a zero value is returned. 
    
    """
    count_list = element_counts(msTuple)
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