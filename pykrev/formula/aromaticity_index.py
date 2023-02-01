from .element_counts import element_counts
import numpy as np
def aromaticity_index(msTuple, index_type = 'rAI'):     
    """ 
    Docstring for function pykrev.aromaticity_index
    ==========
    This function takes an msTuple and returns the aromaticity index of each formula in the formula list.
    
    Use
    ----------
    aromaticity_index(Y)
    
    Returns a numpy.ndarray of len(Y[0]) in which each item, i ,is the aromaticity index corresponding to Y[0][i].  
    
    Parameters
    ----------
    Y: msTuple OR a list of molecular formula strings
    
    index_type: String, one of:
        -'rAI' - reformulated aromaticity index see Mendelez-Perez et al. (2016)
        -'rAImod'- reformulated aromaticity index with modified oxygen coefficient 
        -'AI' - aromaticity index see Koch and Dittmar (2006) 
        -'AImod'- aromaticity index with modified oxygen coefficient see Koch and Dittmar (2006)
    
    
    Info
    ----------
    Reformulated aromaticity index see Mendelez-Perez et al. (2016): 
    "A reformulated aromaticity index equation under consideration for non-aromatic and non-condensed aromatic cyclic carbonyl compounds."
    Aromaticity index see Koch and Dittmar (2006): From mass to structure: an aromaticity index for high-resolution mass data of natural organic matter"
    """ 
    # Tests   
    assert index_type in ['rAI','rAImod','AI','AImod'], 'supply a valid index type, read doc string for more info'
    # Setup
    count_list = element_counts(msTuple)
    AI_array = np.array([])
    # Main
    for count in count_list: 
        try:
            if index_type == 'rAImod':
                AI_counts = (1 + 1/2*((count['C']*2)-count['H']-count['Cl']-count['F']-(count['O']*1)-(count['S']*2)-count['N']-count['P']))/count['C']
            elif index_type == 'rAI':
                AI_counts = (1 + 1/2*((count['C']*2)-count['H']-count['Cl']-count['F']-(count['O']*2)-(count['S']*2)-count['N']-count['P']))/count['C']
            elif index_type == 'AImod':
                AI_counts = (1 + count['C']-(count['Cl']*0.5)-(count['F']*0.5)-(count['H']*0.5)-(count['O']*0.5)-count['S'])/(count['C'] - (count['O']*0.5) - count['N'] - count['S'] - count['P'])
            elif index_type == 'AI':
                AI_counts = (1 + count['C']-count['O']-count['S']-(count['Cl']*0.5)-(count['F']*0.5)-(count['H']*0.5))/(count['C'] - count['O'] - count['N'] - count['S'] - count['P'])
        except ZeroDivisionError:
            AI_counts = np.nan
            print(f"Warning Zero Division Encountered, NaN value returned. {count}")
        AI_array = np.append(AI_array,AI_counts)
    return AI_array
