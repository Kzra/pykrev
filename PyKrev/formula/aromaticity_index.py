import numpy as np
def aromaticity_index(count_list, index_type = 'mod'):
    
          
    """ 
	Docstring for function pyKrev.aromaticity_index
	====================
	This function takes a list of element counts and gives the aromaticity index.
    
	Use
	----
	aromaticity_index(Y)
    
	Returns a numpy array of len(Y) in which each item is the aromaticity index.  
    
	Parameters
	----------
	Y: A list of dictionary items containing atomic counts as produced by element_counts. 
	Must include C,H,O,S,N,P. 
    
	index_type: 'mod' or 'std', which version of the index to use (see info). 
    
	Info
	----------
	Modified aromaticity index see Mendelez-Perez et al. (2016),
	"A reformulated aromaticity index equation under consideration for non-aromatic 
	and non-condensed aromatic cyclic carbonyl compounds."
	rAImod: (1 + 1/2(2C - H - O - 2S - N - P))/C 
	rAIstd: (1 + 1/2(2C - H - 2O - 2S - N - P))/C 
    """    
    
    assert isinstance(count_list,list),'supply a list of counts given by element_counts()'
    assert all(isinstance(i,dict) for i in count_list), 'supply a list of counts given by element_counts()'
    
    AI_array = np.array([])
    warning = 0 
    
    for count in count_list: 
        try:
            if index_type == 'mod':
                AI_counts = (1 + 1/2*((count['C']*2)-count['H']-count['O']-(count['S']*2)-count['N']-count['P']))/count['C']
            elif index_type == 'std':
                AI_counts = (1 + 1/2*((count['C']*2)-count['H']-(count['O']*2)-(count['S']*2)-count['N']-count['P']))/count['C']
        except KeyError:
            raise Exception('Make sure your counts include C,H,O,S,N,P')
        if AI_counts < 0: 
            warning = 1 
            AI_counts = 0
            
        AI_array = np.append(AI_array,AI_counts)
        
    if warning == 1:
        print('Warning: negative ai counts detected and set to zero.')
        
    return AI_array