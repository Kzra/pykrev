from matplotlib import pyplot as plt
from scipy.stats import binned_statistic_2d
from ..formula import element_ratios 
from ..formula import element_counts


def van_krevelen_histogram (formula_list, x_ratio = 'OC', y_ratio ='HC', **kwargs): 
    
    
    """ 
	Docstring for function PyKrev.van_krevelen_histogram
	====================
	This function takes a list of molecular formula strings and plots a van Krevelen histogram. 
    
	Use
	----
	van_krevelen_histogram(Y)
    
	Returns a density index score if bins are provided as integers (see PyKrev.density_index).   
    
	Parameters
	----------
	Y: A list of molecular formula strings.
    x_ratio: element ratio to plot on x axis, given numerator denominator e.g. 'OC'
    y_ratio: element ratio to plot on y axis, given numerator denominator e.g. 'HC'
	**kwargs for pyplot.hist2d() See: https://matplotlib.org/api/_as_gen/matplotlib.pyplot.hist2d.html.
    
    """
    ratio_list = element_ratios(formula_list,ratios=[x_ratio,y_ratio])
        
    if 'bins' not in kwargs: 
        kwargs['bins'] = 20
        xbins = 20
        ybins = 20 
        d_index = density_index(ratio_list,xbins,ybins)
    elif isinstance(kwargs['bins'],int):
        xbins = kwargs['bins']
        ybins = kwargs['bins']
        d_index = density_index(ratio_list,xbins,ybins)
    elif isinstance(kwargs['bins'][0],int) and isinstance(kwargs['bins'][1],int):
        xbins = kwargs['bins'][0]
        ybins = kwargs['bins'][1]
        d_index = density_index(ratio_list,xbins,ybins)
    else: d_index = None

     
    x_axis = []
    y_axis = []
    
    for ratios in ratio_list:
        x_axis.append(ratios['OC'])
        y_axis.append(ratios['HC'])
   
    plt.hist2d(x_axis,y_axis,**kwargs)
    
    plt.xlabel('Atomic ratio of O/C')
    plt.ylabel('Atomic ratio of H/C')
    
    return d_index
        
    
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