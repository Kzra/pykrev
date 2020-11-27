from ..formula import calculate_mass
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import math

def mass_histogram(formula_list,
                   mz_list = [],
                   method = 'monoisotopic',
                   bin_width = [],
                   summary_statistics = False,
                   **kwargs):
    """ 
	Docstring for function PyKrev.mass_histogram
	====================
	This function takes a list of molecular formula strings and plots a histogram of the atomic masses calculated using    method.  
    
	Use
	----
	mass_histogram(Y)
    
    	Returns a tuple containing the mean, median and standard deviation 
    
	Parameters
	----------
	Y: A list of molecular formula strings.
	method: the method used to calculate formula mass. See pk.calculate_mass for more information. 
            one of
            'monoisotopic'
            'average'
            'nominal'
            or
            'mz'- plot a list of mz values provided explicitly by the user
    bin_width: an integer or float specifying the size of each histogram bin
    mz_list: a list or numpy array of mz values, must be provided to use method 'mz'. 
    summary_statistics: boolean, if true print the mean, median and standard deviation on the chart.
    **kwargs: key word arguments to plt.hist
    
    """
    assert method in ['monoisotopic','average','nominal','mz'], 'Provide a valid method. See docstring for info.'
    if method == 'mz':
        assert len(mz_list) > 1, 'provide a list of mz values'
        mass = mz_list
    else: 
        mass = calculate_mass(formula_list, method = method)
    mass_mean = np.mean(mass)
    mass_median = np.median(mass)
    mass_std = np.std(mass)
    if bin_width:
        assert type(bin_width) == int or type(bin_width) == float, 'Provide a scalar value for bin width'
        n = math.ceil((mass.max() - mass.min())/bin_width)
        plt.hist(x=mass, bins = n, **kwargs)
    else:
        plt.hist(x=mass, **kwargs)
    plt.grid(axis='y', alpha=0.75)
    if method == 'mz':
        plt.xlabel('m/z')
    else:
        plt.xlabel(f"{method[0].upper()}{method[1::]} atomic mass")
    plt.ylabel("Counts")
    if summary_statistics:
    ## add to the upper right of the plot
        plt.annotate(f"$\mu={np.round(mass_mean,2)}$\nm = {np.round(mass_median,2)}\n$\sigma={np.round(mass_std,2)}$", xy=(0.75, 0.75), xycoords='axes fraction')
    return mass_mean, mass_median, mass_std