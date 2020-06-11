def missing_formula(*groups):
    
    """ 
	Docstring for function pyKrev.missing_formula
    
	====================
	This function compares n lists of molecular formula and outputs a dictionary containing the missing formula in each list.
    
	Use
	----
	missing_formula(list_1,..,list_n)
    
	Returns a dictionary in which each key corresponds to an input list (labelled "group 1...n") 
	and the corresponding value is a list containing the missing formula in that group.  
    
	Parameters
	----------
	*groups: n lists of molecular formula. Each item in the list should be a formula string. 
    
    """
   
    molecular_formula = dict()
    missing_molecular_formula = dict()
    
    #firstly read in all the groups 
    group_numbers = range(1,len(groups)+1)
    
    for i,g in zip(group_numbers,groups): 
        molecular_formula['group' + str(i)] = g 
    
    #now compare elements of the groups 
    for i in group_numbers:
        #create a shallow copy of the molecular formula dictionary 
        temp_molecular_formula = dict.copy(molecular_formula)
        #create a set out of the unique formulas in the current group 
        current_group = set(molecular_formula['group' + str(i)])
        #remove the current group from the shallow copy
        del temp_molecular_formula['group' + str(i)]
        #create a set out of the unique formulas in the shallow copy
        other_groups = list(temp_molecular_formula.values()) 
        #find the formula that are shared by all other_groups
        common_groups = set(other_groups[0])
        for og in other_groups: 
            common_groups = common_groups & set(og)  # the & operator returns intersecting values of sets
        #find the common formula that aren't present in the current_group 
        missing_formula  = []
        for formula in common_groups:
            if formula in current_group:
                pass
            else:
                missing_formula.append(formula)
        #add it to the unique_molecular_formula dict 
        missing_molecular_formula['group' + str(i)] = missing_formula
    
    
    return missing_molecular_formula
    