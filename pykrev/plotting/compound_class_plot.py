from ..diversity import compound_class
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

def compound_class_plot(formula_list,
                        mass_list = [],
                        method = 'MSCC',
                        **kwargs):
    """ 
	Docstring for function PyKrev.compound_class_plot
	====================
	This function takes a list of molecular formula strings and plots a bar chart of the compound classes present.  
    
	Use
	----
	compound_class_plot(Y)
    
    Returns None
    
	Parameters
	----------
	Y: A list of molecular formula strings.
    mass_list: a list of mz values to pass to pk.compound_class -> required for the MSCC algorithm.
	method: the element to determine the atomic class. One of: C,H,N,O,S or P.
    **kwargs: key word arguments to plt.bar
    
    """
    compoundClass, cclassCounts = compound_class(formula_list,mass_list,method = method)
    labels = []
    values = []
    for c in cclassCounts:      
        labels.append(str(c))
    for v in cclassCounts.values(): 
        values.append(v)
    x_pos = [i for i, _ in enumerate(labels)]
    plt.bar(x_pos, values, **kwargs)
    plt.xticks(x_pos, labels, rotation = 'vertical')
    plt.xlabel("Compound class")
    plt.ylabel("Counts")
    return 