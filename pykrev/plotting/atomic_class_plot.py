from ..formula import element_counts
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

def atomic_class_plot(formula_list,
                      element = 'O',
                      summary_statistics = False,
                      **kwargs):
    """ 
	Docstring for function PyKrev.atomic_class_plot
	====================
	This function takes a list of molecular formula strings and plots a histogram of the atomic classes based on element.  
    
	Use
	----
	atomic_class_plot(Y)
    
    Returns the fig and ax handles and a tuple containing the mean, median and standard deviation 
    
	Parameters
	----------
	Y: A list of molecular formula strings.
	element: the element to determine the atomic class. One of: C,H,N,O,S or P.
    summary_statistics: boolean, if true print the mean, median and standard deviation on the chart.
    **kwargs: key word arguments to plt.hist
    
    """
    
    count_list = element_counts(formula_list)
    countDF = pd.DataFrame(count_list)
    atom = np.array(countDF[element])
    atom_mean = np.mean(atom)
    atom_median = np.median(atom)
    atom_std = np.std(atom)
    plt.hist(x=atom, **kwargs)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel(f"{element} atom class")
    plt.ylabel("Counts")
    fig = plt.gcf()
    ax = plt.gca()
    if summary_statistics:
    ## add to the upper right of the plot
        plt.annotate(f"$\mu={np.round(atom_mean,2)}$\nm = {np.round(atom_median,2)}\n$\sigma={np.round(atom_std,2)}$", xy=(0.75, 0.75), xycoords='axes fraction')
    return fig, ax, (atom_mean, atom_median, atom_std) 