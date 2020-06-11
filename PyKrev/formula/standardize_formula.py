def standardize_formula(qstr,rstr,formula_list):
        
    """ 
	Docstring for function pyKrev.standardize_formula
	====================
	This function takes a list of formula and replaces any instances of qstr with rstr.
	This is an important pre-processing step if isotope names are present in your molecular formula. 
    
	Use
	----
	element_counts(qtr,rstr,Y)
    
	Returns a list of len(Y) in which each item is an atomic formula with qstr replaced by rstr. 
    
	Parameters
	----------
	Y: A list of atomic formula strings. 
    
	qstr: A string with the element name to replace (e.g. '11B')
    
	rstr: A string with the element name to replace with (e.g. 'B')
    
    
    """
    std_formula_list = []
    for formula in formula_list:
        std_formula_list.append(formula.replace(qstr,rstr,1))
    return std_formula_list
