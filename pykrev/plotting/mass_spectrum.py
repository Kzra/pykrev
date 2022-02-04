from numpy.lib.function_base import _rot90_dispatcher, rot90
from ..diversity.normalise_intensity import normalise_intensity
from ..formula.calculate_mass import calculate_mass
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

def mass_spectrum(msTuple, 
                  method = 'monoisotopic',
                  logTransform = False,
                  stepSize = 5,
                  lineColor = 'g',
                  lineWidth = '0.8',
                  invertedAxis = [],
                  invertedAxisLabel = 'Inverted Axis Label',
                  invertedAxisColor = 'b',
                  invertedAxisLineWidth = 0.8,
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
    logTransform: boolean, log transform the y axis
    stepSize: int, number of decimal places to round masses to 
    lineColor: string, color for primary axis
    lineWidth: int, line width for primary axis
    invertedAxis: np.ndarray, array of values to plot on the negative y axis, should all be positive
    invertedAxisLabel: string, label for negative y axis
    invertedAxisColor: string, color for inverted axis
    invertedAxisLineWidth: int, line width for inverted axis
    **kwargs: key word arguments to plt.plot, must not include color or linewidth arguments.
    """
    #Tests
    ## Set default values for color and linewidth unless they have been given in kwargs
    assert 'color' not in kwargs, 'provide colo(u)r as lineColour'
    assert 'lineWidth' not in kwargs, 'provide linewidth as lineWidth'
    assert method in ['monoisotopic','average','nominal','mz'], 'Provide a valid method. See docstring for info.'
    #Setup
    peak_intensities = msTuple.intensity
    mz_list = msTuple.mz
    if logTransform:
        peak_intensities = np.log(peak_intensities)
    if method == 'mz':
        mass = mz_list
    else: 
        mass = calculate_mass(msTuple, method = method)
    if len(invertedAxis) > 0:
        assert len(invertedAxis) == len(peak_intensities), 'inverted data must be the same length as peak intensity array'
        peak_intensities = [x for _,x,_ in sorted(zip(mass,peak_intensities, invertedAxis))]
        mass = [y for y,_,_ in sorted(zip(mass,peak_intensities, invertedAxis))]
        invertedAxis = [z for _,_,z in sorted(zip(mass,peak_intensities, invertedAxis))]
    else:
        peak_intensities = [x for _,x in sorted(zip(mass,peak_intensities))]
        mass = [y for y,_ in sorted(zip(mass,peak_intensities))]
    assert len(peak_intensities) == len(mass)
    ## create plotting mass and plotting peak intensity variables 
    ### step size should be based on the precision of the masses
    ### we can calculate an appropriate step size by finding the minimum difference between adjacent masses
    ### massDiff = np.append(mass,0)[1::] - np.array(mass)
    ### stepSize = str(min(massDiff[::-1]))[::-1].find('.')
    plotting_mass = np.arange(start=min(mass) - 5,stop = max(mass) + 5, step = 10 ** -stepSize)
    plotting_intensity = np.zeros(len(plotting_mass))
    indices = np.nonzero(np.isin(np.around(plotting_mass,stepSize),np.around(mass,stepSize)))[0]
    assert len(indices) == len(mass), f"{len(indices)} {len(mass)} increase stepSize"
    if len(invertedAxis) > 0:
        plotting_inverted = np.zeros(len(plotting_mass))
        for idx,iy,inv in zip(indices,peak_intensities, invertedAxis):
            plotting_intensity[idx] = iy
            # this approach only works with positive values
            if inv < 0:
                inv = 0
            plotting_inverted[idx] = inv
    else:
        for idx,iy in zip(indices,peak_intensities):
            plotting_intensity[idx] = iy
    #Main
    fig, ax1 = plt.subplots()
    ax1.plot(plotting_mass,plotting_intensity,color=lineColor, linewidth=lineWidth, **kwargs)
    if len(invertedAxis) > 0: 
        ax2 = ax1.twinx()
        ax2.plot(plotting_mass,plotting_inverted,color=invertedAxisColor, linewidth=invertedAxisLineWidth, **kwargs)
        ax1.set_ylabel("Intensity", loc = "top")
        ax2.set_ylabel(f"{invertedAxisLabel}", loc = "bottom")
        ax2.grid(axis='y', alpha=0.75)
        ## The following is a bit hacky but I couldn't find a better way. Essentially mirror the axes, and remove the negative tick marks.
        ### First remove negative ticks 
        ax1ticks = [tick for tick in ax1.get_yticks() if tick >=0]
        ax2ticks = [tick for tick in ax2.get_yticks() if tick >=0]
        ### Add a mirrored max negative value
        ax1.set_ylim(ymin=-max(ax1ticks), ymax=max(ax1ticks))
        ax2.set_ylim(ymin=-max(ax2ticks), ymax=max(ax2ticks))
        ### Set tick values
        ax1.set_yticks(ax1ticks)
        ax2.set_yticks(ax2ticks)
        ### Invert ax2
        ax2.invert_yaxis()
        ### Remove the tick label and tick marks for this value
        #yticks = ax1.yaxis.get_major_ticks()
        #yticks[-1].label1.set_visible(False)
        #yticks[-2].label1.set_visible(False)
        #yticks[-1].tick1line.set_visible(False)
        #yticks[-2].tick1line.set_visible(False)
        #yticks = ax2.yaxis.get_major_ticks()
        #yticks[-1].label1.set_visible(False)
        #yticks[-2].label1.set_visible(False)
        #yticks[-1].tick1line.set_visible(False)
        #yticks[-2].tick1line.set_visible(False)
    else:
        ax1.set_ylabel("Intensity")
    ax1.grid(axis='y', alpha=0.75)
    if method == 'mz':
        ax1.set_xlabel('m/z')
    else:
        ax1.set_xlabel(f"{method[0].upper()}{method[1::]} atomic mass")
    if len(invertedAxis) > 0: 
        return fig, ax1, ax2
    else:
        return fig,ax1
