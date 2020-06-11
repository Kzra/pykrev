def aromaticity_index(count_list, index_type = 'mod'):
    
          
    """ 
	Docstring for function pyKrev.aromaticity_index
	====================
	This function takes a list of element counts and gives the aromaticity index.
    
	Use
	----
	aromaticity_index(Y)
    
	Returns a list of len(Y) in which each item is the aromaticity index.  
    
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
    AI_list = []
    
    for count in count_list: 
        try:
            if index_type == 'mod':
                AI_counts = (1 + 1/2*((count['C']*2)-count['H']-count['O']-(count['S']*2)-count['N']-count['P']))/count['C']
            elif index_type == 'std':
                AI_counts = (1 + 1/2*((count['C']*2)-count['H']-(count['O']*2)-(count['S']*2)-count['N']-count['P']))/count['C']
        except KeyError:
            raise Exception('Make sure your counts include C,H,O,S,N,P')
        AI_list.append(AI_counts)
    
    return AI_list