from matplotlib import pyplot as plt
import numpy as np
def spiral_plot(msTuple, colour=[], size=[], radius=1, theta=3000, mass_order = 'ascending', colourmap = 'viridis'):
    """ 
    Docstring for function pykrev.spiral_plot
    ==========
    Takes an msTuple as input and creates a spiral molecular formulae plot
    
    Use 
    ----------
    spiral_plot(Y) 
    
    Returns the figure and axis handle for the spiral_plot
    
    Parameters 
    ----------
    Y: An msTuple 
    
    colour: array or string, colour of the points in the scatter plot (see matplotlib 'c' kwarg)
        default is colour by intensity
    size: array or float, size of the points in the scatter plot (see matplotlib 'size' kwarg)
        default is colour by relative intensity
    radius: Int, radius of the spiral
    theta: Int, determines the shape of the spiral (<1000 truncated spiral, 1000-9999 full spiral, >9999 funky spiral)
    mass_order: String, which order to plot the formulae. One of:
        'ascending': lowest masses in centre of spiral, increasing outward
        'descending': lowest masses in outside of spiral, increasing inward
    colourmap: String, cmap used to colour points (see matplotlib 'cmap' kwarg)
    """
    #Tests
    assert(mass_order in ['ascending','descending']), "incorrect mass_order given"
    #Setup
    if len(colour) == 0:
        colour = msTuple.intensity
    if len(size) == 0:
        size = msTuple.intensity/max(msTuple.intensity) * 100
    N = len(msTuple.formula)
    if mass_order == 'ascending':
        sortIdx = np.argsort(msTuple.mz)
    elif mass_order == 'descending':
        sortIdx = np.argsort(msTuple.mz)[::-1]
    colour = np.take_along_axis(colour, sortIdx, axis=0)
    size = np.take_along_axis(size, sortIdx, axis=0)
    r = np.linspace(0,radius,N)
    t = np.linspace(0,theta,N)
    x = r*np.cos(np.radians(t))
    y = r*np.sin(np.radians(t))
    #Main
    fig, ax = plt.subplots()
    ax.scatter(x,y,c=colour,s=size, cmap=colourmap)
    ax.set_axis_off()
    return fig, ax
    