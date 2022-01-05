import itertools
import numpy as np
import pandas as pd
def find_intersections(msTupleDict, exclusive = True):
    """   
	Docstring for function pykrev.find_intersections
    
	====================
	This function compares n msTuples inside an msTupleDict and outputs a dictionary containing the intersections between each msTuple.
    
	Use
	----
	find_intersections(Y)
    
	Returns a dictionary in which each key corresponds to a combination of sample names 
	and the corresponding value is a set containing the intersections between the groups in that combination.  
    
	Parameters
	----------
	Y: msTupleDict
	exclusive: Boolean, True or False, depending on whether you want the intersections to contain only unique values.
    """
    #Setup
    formula_lists = []
    group_labels = msTupleDict.keys()
    for msTuple in msTupleDict.values():
        formula_lists.append(msTuple.formula)
    combinations = [seq for i in range(0,len(group_labels)+1) for seq in itertools.combinations(group_labels,i) if len(seq) > 0]
    combinations = sorted(combinations,key = lambda c : len(c),reverse = True) # sort combinations by length
    if exclusive == True:
        assigned_formula = set() #create a set that will hold all the formula already assigned to a group
    amb = pd.DataFrame(data = formula_lists).T
    amb.columns = group_labels
    intersections = dict() 
    #Main
    for combo in combinations:         
            queries = [] 
            for c in combo:
                formula = list(filter(None,amb[c])) #Remove None entries introduced by dataframe
                queries.append(set(formula)) 
            if len(queries) == 1: #if there is only one query find the unique elements in it 
                q_set = frozenset(queries[0]) #qset is a frozen set, so it will not be mutated by changes to queries[0]
                for f_list in formula_lists: #cycle all formula in formula_lists
                    set_f = frozenset(f_list) #convert f_list to sets, must be frozen so type matches q_set
                    if set_f == q_set: # ignore the set that corresponds to the query 
                        pass
                    else:
                        queries[0] = queries[0] - set_f #delete any repeated elements in fset
                            
                intersections[combo] = queries[0] 
                
            elif len(queries) > 1: 
                if exclusive == True:
                    q_intersect = intersect(queries)
                    intersections[combo] = q_intersect - assigned_formula #remove any elements from q_intersect that have already been assigned
                    assigned_formula.update(q_intersect) #update the assigned_set with q_intersect
                else:
                    intersections[combo] = intersect(queries)  
    return intersections
    
def intersect(samples,counter=0):
    """ This command uses recursion to find the intersections between a variable number of sets given in samples. 
        Where samples = [set_1,set_2,...,set_n] """
    if len(samples) == 1: 
        return samples[0]
    a = samples[counter]
    b = samples[counter+1::]
    if len(b) == 1: #check to see whether the recursion has reached the final element
        return a & b[0]
    else:
        counter += 1 
        return a & intersect(samples,counter)