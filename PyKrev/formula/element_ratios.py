def element_ratios(count_list, ratios = ['OC','HC']):
        
    
    """ 
	Docstring for function pyKrev.element_ratios
	====================
	This function takes a list of element counts and gives atomic ratios for each ratio given in ratios.
    
	Use
	----
	element_ratios(Y)
    
	Returns a list of len(Y) in which each item is a dictionary containing the atomic ratios listed in ratios, 
	calculated as the first element / second element. 
    
	Parameters
	----------
	Y: A list of dictionary items containing atomic counts as produced by element_counts.
    
	ratios: A list of atomic ratios to compute, calculated as first element divided by second element. 
	Only C,H,N,O,P & S can be used. 
    
        """            
    ratio_list = []
    
    for count in count_list:
        ratio_counts = dict()
        
        for atom in ratios:
            if len(atom) > 2:
                raise Exception('Only C,H,N,O,P & S can be used. ')
            
            ratio_counts[atom] = count[atom[0]]/count[atom[1]]
        ratio_list.append(ratio_counts)
        
    return ratio_list