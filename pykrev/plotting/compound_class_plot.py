from ..diversity.compound_class import compound_class
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

def compound_class_plot(msTuple,
                        method = 'MSCC',
                        **kwargs):
    """ 
	Docstring for function PyKrev.compound_class_plot
	====================
	This function takes an mStuple and plots a bar chart of the compound classes present.  
    
	Use
	----
	compound_class_plot(Y)
    
    Returns the figure and axes handles and a tuple containing a list of compound class assignments and a dictionary containing the compound class counts
    
	Parameters
	----------
	Y: msTuple
    mass_list: a list of mz values to pass to pk.compound_class -> required for the MSCC algorithm.
	method: the element to determine the atomic class. One of: C,H,N,O,S or P.
    **kwargs: key word arguments to plt.bar
    """
    compoundClass, cclassCounts = compound_class(msTuple,method = method)
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
    fig = plt.gcf()
    ax = plt.gca()
    return fig, ax, (compoundClass, cclassCounts)