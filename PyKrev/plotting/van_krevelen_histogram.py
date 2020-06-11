from matplotlib import pyplot as plt
from scipy.stats import binned_statistic_2d
from matplotlib import cm

def van_krevelen_histogram (ratio_list,colourmap = 'CMRmap',xbins = 20, ybins = 20): 
    
    
    """ 
	Docstring for function PyKrev.van_krevelen_histogram
	====================
	This function takes a list of H/C-O/C ratios and plots a van Krevelen histogram. 
    
	Use
	----
	van_krevelen_histogram(Y)
    
	Returns figure, axes, and a density index score (see PyKrev.density_index).   
    
	Parameters
	----------
	Y: A list of atom ratios (must contain H/C and O/C). See PyKrev.element_ratios.
    
	colourmap: colours to map the density values to. See: https://matplotlib.org/3.2.1/tutorials/colors/colormaps.html
    
	xbins: number of x axis bins. See: https://matplotlib.org/api/_as_gen/matplotlib.pyplot.hist2d.html
	ybins: number of y axis bins. 
    
    """
    
    figure = plt.figure()
    
    x_axis = []
    y_axis = []
    
    for ratios in ratio_list:
        x_axis.append(ratios['OC'])
        y_axis.append(ratios['HC'])
   
    histogram = plt.hist2d(x_axis,y_axis,(xbins,ybins), cmap = colourmap)
    d_index = density_index(ratio_list,xbins,ybins)
    
    cbar = plt.colorbar()
    cbar.set_label('Number of moelcular formula')
    
    plt.xlabel('Atomic ratio of O/C')
    plt.ylabel('Atomic ratio of H/C')
    
    return histogram, d_index
        
    
def density_index (ratio_list,xbins=20,ybins=20):
    
    """      
	Docstring for function PyKrev.van_krevelen_histogram
	====================
	This function takes a list of H/C-O/C ratios and calculates a density score. 
	This is achieved by dividing the average number of points by the number of points in the most populated. 
	Giving average relative density. For 100 bins a score of 1 means all bins are equally dispersed, and a score of 0.01 
	means all points fall into one bin.
    
	Use
	----
	density_index(Y)
    
	Returns a float containing density index score.   
    
	Parameters
	----------
	Y: A list of atom ratios (must contain H/C and O/C). See PyKrev.element_ratios.
    
	xbins: number of x axis bins. See: https://matplotlib.org/api/_as_gen/matplotlib.pyplot.hist2d.html
	ybins: number of y axis bins. 
    
    """
    
    x = []
    y = []
    
    for ratios in ratio_list:
        x.append(ratios['OC'])
        y.append(ratios['HC'])
        
    #count the number of points that fall into each of our pre-defined bins
    bin_results = binned_statistic_2d(x,y,None,'count',bins = (xbins,ybins))
    bin_counts = bin_results.statistic
    #compute the average bin density
    average_bin_density = bin_counts.mean()
    #compute the max bin density 
    max_bin_density = bin_counts.max()
    #compute the density distribution 
    d_index = average_bin_density/max_bin_density
          
    return d_index