def element_counts(msTuple):
    """ 
	Docstring for function pykrev.element_counts
	====================
	This function takes an msTuple and gives atomic counts for C,H,N,O,P,S,Cl & F for each formula in the formula list.
    
	Use
	----
	element_counts(Y)
    
	Returns a list of len(Y) in which each element, i , is a dictionary containing the atomic counts of Y[0][i]. 
    
	Parameters
	----------
	Y: msTuple OR a list of elemental formula strings. 
    
    Info
    ----------
    For formula strings:
        All integers must be standard script (e.g. C6H8O7). 
        Formula strings should not contain isotopologues (e.g. C9H12O6 13C1) 
        Formula strings should only contain C,H,N,O,P,S,Cl and F atoms.
    """
    #Tests
    if type(msTuple) == list:
        formula_list = msTuple
    else:
        formula_list = msTuple[0]
    #Setup
    elements=['C','H','N','O','P','S','Cl','F']
    count_list = []
    #Main
    for formula in formula_list:
        element_numbers = dict()
        for element in elements:
            element_numbers[element] = 0
            alpha_idx = formula.find(element) #find this first instance of the element name. If the formula is given with alphabetically similar two digit element names e.g. Co before single digit names e.g. C. 
                                           #it won't work. In theory this shouldn't happen because formulas should be listed alphabetically, therefore single character names will come before two character names. 
            if alpha_idx > -1: # if the element was found in the formula   
                element_numbers[element] = 1 #it must have at least 1 atom 
                for i in range(alpha_idx+len(element)+1,(len(formula)+1)):
                    try:
                        element_numbers[element] = int(formula[alpha_idx+len(element):i]) #Count the number of atoms belonging to this element starting with the next character, and progressively making the character window bigger
                    except ValueError: # this occurs when you have reached the next element name
                        break
        count_list.append(element_numbers)
    return count_list
