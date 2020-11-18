import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import gaussian_kde
from ..formula import element_ratios 

def van_krevelen_plot(formula_list,
                      x_ratio = 'OC',
                      y_ratio = 'HC',
                      **kwargs):
    
    """ 
	Docstring for function PyKrev.van_krevelen_plot
	====================
	This function takes a list of molecular formula strings and plots a van Krevelen diagram. 
    
	Use
	----
	van_krevelen_plot(Y)
  
    
	Parameters
	----------
	Y: A list of molecular formula strings.
	colour: A list or numpy array of floats or integers of len(Y) or
        'density' : kernel density see https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html
    x_ratio: element ratio to plot on x axis, given numerator denominator e.g. 'OC'
    y_ratio: element ratio to plot on y axis, given numerator denominator e.g. 'HC'
    **kwargs: key word arguments for pyplot.scatter(). 

	Info
	----------
	Given a list of ratios this function will produce Van Krevelen plots as described in 
	Kim, Sunghwan, Robert W. Kramer, and Patrick G. Hatcher. 
	"Graphical method for analysis of ultrahigh-resolution broadband mass spectra of natural organic matter, 
	the van Krevelen diagram." 
	 Analytical Chemistry 75.20 (2003): 5336-5344. 

    """
    ratio_list = element_ratios(formula_list,ratios=[x_ratio,y_ratio])
    #check that color is provided as 'c'
    assert 'color' not in kwargs, 'supply key word color as c'
        
    if 'c' not in kwargs: 
        kwargs['c'] = ['blue'] * len(ratio_list)
    elif isinstance(kwargs['c'],str) and kwargs['c'] == 'density':
        kwargs['c'] = kernel_density(ratio_list)
    elif len(kwargs['c']) != len(ratio_list):
        raise ValueError('colour list and ratio list must be the same length.')

    x_axis = []
    y_axis = []

    for ratios in ratio_list:
            x_axis.append(ratios[x_ratio])
            y_axis.append(ratios[y_ratio])
            
            
    plt.scatter(x_axis, y_axis, **kwargs)
    
    #apply grid lines 
    plt.grid(True) 
        
    #add on chemical reaction lines 
    #slopes are taken from Hatcher et al. (2003) Graphical method for analysis...
    #if 'hydrogenation' in ref_reactions: 
        #plt.plot((0.5,0.5),(2.0,0.5), "r--",alpha=1)
        #plt.text(0.52,0.52,'H',fontsize=10,alpha=1,color='r')
    #if 'redox' in ref_reactions:
        #plt.plot((0.1,0.8),(1,1), "r--",alpha=1)
        #plt.text(0.81,1.02,'R/O',fontsize=10,alpha=1,color='r')
    #if 'condensation' in ref_reactions:
        #plt.plot((0.2,0.8),(0.4,1.6), "r--",alpha=1)
        #plt.text(0.82,1.58,'Co',fontsize=10,alpha=1,color='r')
    #if 'methylation' in ref_reactions: 
        #plt.plot((0.1,.6),(1.8,.8),"r--",alpha = 1)
        #plt.text(.61,.79,'Me',fontsize=10,alpha= 1,color='r')
    #if 'carboxylation' in ref_reactions:
        #plt.plot((0.1,0.8),(2,2),"r--",alpha=1)
        #plt.text(0.82,2.02,'Cbx',fontsize=10,alpha=1,color='r')
        

    plt.xlabel(f"Atomic ratio of {x_ratio[0]}/{x_ratio[1]}")
    plt.ylabel(f"Atomic ratio of {y_ratio[0]}/{y_ratio[1]}")
    
    return 

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