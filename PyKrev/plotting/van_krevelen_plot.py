import numpy as np

from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

from scipy.stats import gaussian_kde

def van_krevelen_plot(ratio_list,colour=[],scale=[], symbol='o',colourmap ='viridis', max_marker_size = 20, transparency = 0.5,
                        ref_compounds = ['lipid','condensed hydrocarbon', 'lignin', 'protein','cellulose'], 
                        ref_reactions = ['methylation','hydrogenation','condensation','redox','carboxylation']):
    
    """ 
	Docstring for function PyKrev.van_krevelen_plot
	====================
	This function takes a list of H/C-O/C ratios and plots a van Krevelen diagram. 
    
	Use
	----
	van_krevelen_plot(Y)
    
	Returns figure, axes, colour bar and legend handles for a van krevelen plot.   
    
	Parameters
	----------
	Y: A list of atom ratios (must contain H/C and O/C). See PyKrev.element_ratios.
    
	colour: A list or numpy array of floats or integers of len(Y) or
	 'density' : kernel density see https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html

	scale: a list or numpy array of floats or integers of len(Y)
    
	symbol: scatter plot symbol to use. See: https://matplotlib.org/3.2.1/api/markers_api.html
    
	colourmap: colours to map the values in 'colour' too. See: https://matplotlib.org/3.2.1/tutorials/colors/colormaps.html
    
	max_marker_size: maximum size of the values in 'scale'.
    
	trasparency: transparency of the marker points. 
    
	ref_reactions: a list of strings detailing reaction lines to include on the van krevelen plot - 
	 'methylation' : Me
	 'hydrogenation' : H
	 'condenstation' : Co
	 'redox' : R/O
	 'carboxylation' : Cbx    
     
	ref_compounds: a list of strings  

	Info
	----------
	Given a list of ratios this function will produce Van Krevelen plots as described in 
	Kim, Sunghwan, Robert W. Kramer, and Patrick G. Hatcher. 
	"Graphical method for analysis of ultrahigh-resolution broadband mass spectra of natural organic matter, 
	the van Krevelen diagram." 
	 Analytical Chemistry 75.20 (2003): 5336-5344. 
         
    
    """ 
        
    x_axis = []
    y_axis = []
    cbar = []
    legend = []
    
    fig, ax = plt.subplots()
    
    for ratios in ratio_list:
        x_axis.append(ratios['OC'])
        y_axis.append(ratios['HC'])
        
    if not scale:
        std_scale = [max_marker_size]*len(ratio_list)
    
    else:
        
        std_scale = np.zeros(len(scale))
        max_scale = np.max(scale)
        
        for i,n in zip(range(0,len(std_scale)),scale):
            std_scale[i] = n/max_scale * max_marker_size #standardize scale values 
            
    
    if colour == 'density':
        colour = kernel_density(ratio_list)
    
    if not colour: 
        std_color = ['blue'] * len(ratio_lists)
        
    else:
        
        std_color = np.zeros(len(colour))
        max_color = np.max(colour)
        
        for i,n in zip(range(0,len(std_color)),colour):
            std_color[i] = n/max_color  #standardize colour valeues
    
        
    scatter = ax.scatter(x_axis, y_axis, alpha=transparency, edgecolors='none', s=std_scale, c=std_color, marker = symbol, cmap = colourmap)
    
    #apply grid lines 
    ax.grid(True) 
    
    # produce a legend with a cross section of sizes from the scatter
    if scale:
        kw = dict(prop="sizes",func=lambda std_scale: std_scale/max_marker_size * max_scale,
                  color = scatter.cmap(0.5),num = 5, alpha = 0.6)  

        legend = ax.legend(*scatter.legend_elements(**kw),
                        loc="lower right")
    
    #construct a colorbar, normalise the values in colour 
    if colour: #user has supplied a list of colours
        scalarmappable = cm.ScalarMappable(cmap=colourmap) 
        scalarmappable.set_array(colour)
        cbar = plt.colorbar(scalarmappable)
    
    #add on chemical reaction lines 
    #slopes are taken from Hatcher et al. (2003) Graphical method for analysis...
    if 'hydrogenation' in ref_reactions: 
        ax.plot((0.5,0.5),(2.0,0.5), "r--",alpha=0.7)
        ax.text(0.52,0.52,'H',fontsize=10,alpha=0.7,color='r')
    if 'redox' in ref_reactions:
        ax.plot((0.1,0.8),(1,1), "r--",alpha=0.7)
        ax.text(0.81,1.02,'R/O',fontsize=10,alpha=0.7,color='r')
    if 'condensation' in ref_reactions:
        ax.plot((0.2,0.8),(0.4,1.6), "r--",alpha=0.7)
        ax.text(0.82,1.58,'Co',fontsize=10,alpha=0.7,color='r')
    if 'methylation' in ref_reactions: 
        ax.plot((0.1,.6),(1.8,.8),"r--",alpha = 0.7)
        ax.text(.61,.79,'Me',fontsize=10,alpha=0.57,color='r')
    if 'carboxylation' in ref_reactions:
        ax.plot((0.1,0.8),(2,2),"r--",alpha=0.7)
        ax.text(0.82,2.02,'Cbx',fontsize=10,alpha=0.7,color='r')
        
    #add on compound class polygons 
    ##compound classes are taken from ...
    

    plt.xlabel('Atomic ratio of O/C')
    plt.ylabel('Atomic ratio of H/C')
    
    
    return fig, ax, cbar, legend 

def kernel_density(ratio_list):
    
    """This function computes the kernel density of a list of molecular formula 'OC' and 'HC ratios using gaussian kernels.
       It returns a list containing the corresponding density values. For information on this function see 
       https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html                        """ 
        
    x = []
    y = []
    
    for ratios in ratio_list:
        x.append(ratios['OC'])
        y.append(ratios['HC'])
    
    xy =  np.vstack([x,y])
    kd = gaussian_kde(xy)(xy) #calling the inner function on (xy) and then the result on (xy)
    
    return list(kd)