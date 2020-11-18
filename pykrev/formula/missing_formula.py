def missing_formula(*groups,group_labels = []):
    
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
	 group_labels = list of group label strings of len(groups)
    
    """
    
    default_names = ['group 1','group 2','group 3','group 4','group 5','group 6','group 7','group 8','group 9']

    if not group_labels:
        group_labels = default_names[0:len(groups)]
   
    molecular_formula = dict()
    missing_molecular_formula = dict()
    
    for i,g in zip(group_labels,groups): 
        molecular_formula[i] = g 
    
    #now compare elements of the groups 
    for i in group_labels:
        #create a shallow copy of the molecular formula dictionary 
        temp_molecular_formula = dict.copy(molecular_formula)
        #create a set out of the unique formulas in the current group 
        current_group = set(molecular_formula[i])
        #remove the current group from the shallow copy
        del temp_molecular_formula[i]
        #create a list out of the formulas in the shallow copy
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
        missing_molecular_formula[i] = missing_formula
    
    
    return missing_molecular_formula
    