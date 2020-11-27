from ..diversity import relative_intensity
from ..formula import calculate_mass
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

def mass_spectrum(formula_list, 
                  peak_intensities,
                  mz_list = [],
                  normalise = True,
                  method = 'monoisotopic',
                  **kwargs):
    """ 
	Docstring for function PyKrev.mass_spectrum
	====================
	This function takes a list of molecular formula strings and a corresponding numpy array of peak intensities and plots a     mass spectrum using atomic masses calculated via method.  
    
	Use
	----
	mass_spectrum(Y,peak_intensities)
    
    	Returns None
    
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
    peak_intensities: a list or numpy array of peak intensities 
    normalise: boolean, whether to calculate relative intensity values 
    mz_list: a list or numpy array of mz values, must be provided to use method 'mz'.
    summary_statistics: boolean, if true print the mean, median and standard deviation on the chart.
    **kwargs: key word arguments to plt.plot
    
    """
    ## Set default values for color and linewidth unless they have been given in kwargs
    assert 'c' not in kwargs, 'provide colo(u)r as color'
    if 'color' not in kwargs:
        color = 'black'
    else:
        color = kwargs['color']
    if 'linewidth' not in kwargs:
        linewidth = 0.8
    else:
        linewidth = kwargs['linewidth'] 
    assert len(peak_intensities) > 0, 'Must provide a list or numpy array of peak intensities'
    assert len(peak_intensities) == len(formula_list), 'peak_intensities and formula list must be the same length'
    assert method in ['monoisotopic','average','nominal','mz'], 'Provide a valid method. See docstring for info.'
    if normalise:
        peak_intensities = relative_intensity(peak_intensities)
    if method == 'mz':
        assert len(mz_list) > 1, 'provide a list of mz values'
        mass = mz_list
    else: 
        mass = calculate_mass(formula_list, method = method)
    # perform sorting of the peak_intensity and mass variables to ensure they are in ascending mass order
    peak_intensities = [x for _,x in sorted(zip(mass,peak_intensities))]
    mass = [y for y,_ in sorted(zip(mass,peak_intensities))]
    plt.plot(mass,peak_intensities,color=color, linewidth=linewidth)
    plt.grid(axis='y', alpha=0.75)
    if method == 'mz':
        plt.xlabel('m/z')
    else:
        plt.xlabel(f"{method[0].upper()}{method[1::]} atomic mass")
    plt.ylabel("Intensity")
    return 