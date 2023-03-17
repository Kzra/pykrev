import numpy as np
from .element_counts import element_counts
def element_ratios(msTuple, ratios = ['OC','HC']):
        
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
    
	ratios: list, atomic ratios to compute, calculated as first element divided by second element. 
        
    """
    #Tests
    for ratio in ratios:
        assert len(ratio) < 4 and len(ratio) > 1, 'Too many or too few characters in the ratio'
    #Setup
    count_list = element_counts(msTuple)
    ratio_list = []
    #Main
    for count in count_list:
        ratio_counts = dict()
        for ratio in ratios:
            # we need to set the second index to 2 if chlorine is the first element in the ratio  as it has two characters in the element name
            if 'Cl' in ratio and ratio.index('Cl') == 0:
                sidx = 2
            else:
                sidx = 1
            # this if statement checks for division by zero and sets value to nan if so
            if count[ratio[sidx]] == 0:
                    ratio_counts[ratio] = np.nan
            else:
                ratio_counts[ratio] = count[ratio[0:sidx]]/count[ratio[sidx::]]
        ratio_list.append(ratio_counts)    
    return ratio_list