from .element_counts import element_counts
import numpy as np
def aromaticity_index(formula_list, index_type = 'rAI'):
    
          
    """ 
	Docstring for function pyKrev.aromaticity_index
	====================
	This function takes a list of molecular formula strings and returns the aromaticity index.
    
	Use
	----
	aromaticity_index(Y)
    
	Returns a numpy array of len(Y) in which each item is the aromaticity index.  
    
	Parameters
	----------
	Y: A list of molecular formula strings. 
    
	index_type: One of the following strings:
        -'rAI' - reformulated aromaticity index see Mendelez-Perez et al. (2016)
        -'rAImod'- reformulated aromaticity index with modified oxygen coefficient 
        -'AI' - aromaticity index see Koch and Dittmar (2006) 
        -'AImod'- aromaticity index with modified oxygen coefficient see Koch and Dittmar (2006)
    
	Info
	----------
	Reformulated aromaticity index see Mendelez-Perez et al. (2016),
	"A reformulated aromaticity index equation under consideration for non-aromatic 
	and non-condensed aromatic cyclic carbonyl compounds."
    Aromaticity index see Koch and Dittmar (2006)
    "From mass to structure: an aromaticity index for
     high-resolution mass data of natural organic matter"

    """    
    
    count_list = element_counts(formula_list)
    assert index_type in ['rAI','rAImod','AI','AImod'], 'supply a valid index type, read doc string for more info'
    AI_array = np.array([])
    warning = 0 
    
    for count in count_list: 
        try:
            if index_type == 'rAImod':
                AI_counts = (1 + 1/2*((count['C']*2)-count['H']-(count['O']*1)-(count['S']*2)-count['N']-count['P']))/count['C']
            elif index_type == 'rAI':
                AI_counts = (1 + 1/2*((count['C']*2)-count['H']-(count['O']*2)-(count['S']*2)-count['N']-count['P']))/count['C']
            elif index_type == 'AImod':
                AI_counts = (1 + count['C'] - (count['H']*0.5)-(count['O']*0.5)- count['S'])/(count['C'] - (count['O']*0.5) - count['N'] - count['S'] - count['P'])
            elif index_type == 'AI':
                AI_counts = (1 + count['C'] -count['O']- count['S'] - (0.5*count['H']))/(count['C'] - count['O'] - count['N'] - count['S'] - count['P'])
        except ZeroDivisionError:
            AI_counts = np.nan
            print(f"Warning Zero Division Encountered, NaN value returned. {count}")
        if AI_counts < 0: 
            warning = 1 
            AI_counts = 0
            
        AI_array = np.append(AI_array,AI_counts)
        
    if warning == 1:
        print('Warning: negative ai counts detected and set to zero.')
        
    return AI_array