from ..diversity.normalise_intensity import normalise_intensity
from ..formula.calculate_mass import calculate_mass
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

def mass_spectrum(msTuple, 
                  normalise = True,
                  method = 'monoisotopic',
                  **kwargs):
    """ 
	Docstring for function PyKrev.mass_spectrum
	====================
	This function takes an msTuple and plots a mass spectrum using atomic masses calculated via method.  
    
	Use
	----
	mass_spectrum(Y,peak_intensities)
    
    Returns the figure and axes handles
    
	Parameters
	----------
	Y: msTuple
	method: the method used to calculate formula mass. See pk.calculate_mass for more information. 
            one of
            'monoisotopic'
            'average'
            'nominal'
            or
            'mz'- plot a list of mz values in msTuple
    normalise: boolean, whether to calculate relative intensity values 
    summary_statistics: boolean, if true print the mean, median and standard deviation on the chart.
    **kwargs: key word arguments to plt.plot
    """
    #Tests
    ## Set default values for color and linewidth unless they have been given in kwargs
    assert 'c' not in kwargs, 'provide colo(u)r as color'
    assert method in ['monoisotopic','average','nominal','mz'], 'Provide a valid method. See docstring for info.'
    #Setup
    if 'color' not in kwargs:
        color = 'black'
    else:
        color = kwargs['color']
    if 'linewidth' not in kwargs:
        linewidth = 0.8
    else:
        linewidth = kwargs['linewidth'] 
    peak_intensities = msTuple[1]
    mz_list = msTuple[2]
    if normalise:
        peak_intensities = normalise_intensity(peak_intensities)
    if method == 'mz':
        mass = mz_list
    else: 
        mass = calculate_mass(msTuple, method = method)
    peak_intensities = [x for _,x in sorted(zip(mass,peak_intensities))]
    mass = [y for y,_ in sorted(zip(mass,peak_intensities))]
    #Main
    plt.plot(mass,peak_intensities,color=color, linewidth=linewidth)
    plt.grid(axis='y', alpha=0.75)
    if method == 'mz':
        plt.xlabel('m/z')
    else:
        plt.xlabel(f"{method[0].upper()}{method[1::]} atomic mass")
    plt.ylabel("Intensity")
    fig = plt.gcf()
    ax = plt.gca()
    return fig,ax
