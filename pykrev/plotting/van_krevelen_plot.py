import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import gaussian_kde
from ..formula.element_ratios import element_ratios 
from matplotlib.patches import Rectangle


def van_krevelen_plot(msTuple,
                      x_ratio = 'OC',
                      y_ratio = 'HC',
                      patch_classes = [],
                      patch_alpha = 0.3,
                      patch_colors = ['#762a83','#9970ab','#c2a5cf','#e7d4e8','#d9f0d3','#a6dba0','#5aae61','#1b7837'],
                      patch_text = True,
                      **kwargs):
    
    """ 
	Docstring for function PyKrev.van_krevelen_plot
	====================
	This function takes an msTuple and plots a van Krevelen diagram. 
    
	Use
	----
	van_krevelen_plot(Y)
  
    Returns the figure and axes handles of the plot 
    
	Parameters
	----------
	Y: msTuple
	colour: A list or numpy array of floats or integers of len(Y[0]) or
        'density' : kernel density see https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html
    x_ratio: element ratio to plot on x axis, given numerator denominator e.g. 'OC'
    y_ratio: element ratio to plot on y axis, given numerator denominator e.g. 'HC'
    patch_classes: a list of the compound classes boundaries (taken from formularity software) to overlay as patches, can include:
        'lipid-like'
        'carbohydrate-like'
        'unsaturated hydrocarbons'
        'condensed aromatics'
        'lignin-like'
        'tannin-like'
        'amino sugar-like'
        'protein-like'
    patch_alpha: the transparency of the compound class  (float between 0 and 1)
    patch_colors: hex values for the colors of each class in patch_classes 
    patch_text: boolean, include text labels on the compound class patches 
    **kwargs: key word arguments for pyplot.scatter(). 

	Info
	----------
	Given a list of ratios this function will produce Van Krevelen plots as described in 
	Kim, Sunghwan, Robert W. Kramer, and Patrick G. Hatcher. 
	"Graphical method for analysis of ultrahigh-resolution broadband mass spectra of natural organic matter, 
	the van Krevelen diagram."  
	Analytical Chemistry 75.20 (2003): 5336-5344. 
    """
    #Tests
    #check that color is provided as 'c'
    assert 'color' not in kwargs, 'supply key word color as c'
    ratio_list = element_ratios(msTuple,ratios=[x_ratio,y_ratio])
    if 'c' not in kwargs: 
        kwargs['c'] = ['blue'] * len(ratio_list)
    elif isinstance(kwargs['c'],str) and kwargs['c'] == 'density':
        kwargs['c'] = kernel_density(ratio_list)
    elif len(kwargs['c']) != len(ratio_list):
        raise ValueError('colour list and ratio list must be the same length.')
    assert len(patch_colors) >= len(patch_classes), "Provide at least as many colors as classes"
    #Setup
    x_axis = []
    y_axis = []
    for ratios in ratio_list:
            x_axis.append(ratios[x_ratio])
            y_axis.append(ratios[y_ratio])
    #Main
    plt.scatter(x_axis, y_axis, **kwargs)
    #apply grid lines 
    plt.grid(True) 
    #add on chemical class patches
    #boundaries taken from formularity software
    cindx = 0 #index for the patch colours 
    if 'lipid-like' in patch_classes:
        plt.gca().add_patch(Rectangle((0.01,1.5),0.29,0.7,linewidth=2,edgecolor =patch_colors[cindx],facecolor=patch_colors[cindx],alpha = patch_alpha))
        if patch_text: plt.text(0.012,2.14,'Lipid',fontsize=11,alpha=1, color = 'k')
        cindx += 1 
    if 'carbohydrate-like' in patch_classes:
        plt.gca().add_patch(Rectangle((0.7,1.5),0.4,0.8,linewidth=2,edgecolor =patch_colors[cindx],facecolor=patch_colors[cindx],alpha = patch_alpha))
        if patch_text: plt.text(0.702,2.24,'Carbs',fontsize=11,alpha=1, color = 'k')
        cindx += 1 
    if 'unsaturated hydrocarbons' in patch_classes:
        plt.gca().add_patch(Rectangle((0.01,0.8),0.09,0.7,linewidth=2,edgecolor =patch_colors[cindx],facecolor=patch_colors[cindx],alpha = patch_alpha))
        if patch_text: plt.text(0.012,1.44,'Unsat HC',fontsize=11,alpha=1, color = 'k')
        cindx += 1
    if 'condensed aromatics' in patch_classes:
        plt.gca().add_patch(Rectangle((0.01,0.2),0.09,0.6,linewidth=2,edgecolor =patch_colors[cindx],facecolor=patch_colors[cindx],alpha = patch_alpha))
        if patch_text: plt.text(0.012,0.74,'Con HC',fontsize=11,alpha=1, color = 'k')
        cindx += 1 
    if 'lignin-like' in patch_classes:
        plt.gca().add_patch(Rectangle((0.1,0.8),0.6,0.8,linewidth=2,edgecolor =patch_colors[cindx],facecolor=patch_colors[cindx],alpha = patch_alpha))
        if patch_text: plt.text(0.102,1.54,'Lignin',fontsize=11,alpha=1, color = 'k')
        cindx += 1
    if 'tannin-like' in patch_classes:
        plt.gca().add_patch(Rectangle((0.7,0.8),0.5,0.8,linewidth=2,edgecolor =patch_colors[cindx],facecolor=patch_colors[cindx],alpha = patch_alpha))
        if patch_text: plt.text(0.702,1.54,'Tannin',fontsize=11,alpha=1, color = 'k')
        cindx += 1
    if 'amino sugar-like' in patch_classes:
        plt.gca().add_patch(Rectangle((0.6,1.5),0.1,0.7,linewidth=2,edgecolor =patch_colors[cindx],facecolor=patch_colors[cindx],alpha = patch_alpha))
        if patch_text: plt.text(0.602,2.14,'AminoSugar',fontsize=11,alpha=1, color = 'k')
        cindx += 1 
    if 'protein-like' in patch_classes:
        plt.gca().add_patch(Rectangle((0.3,1.5),0.3,0.8,linewidth=2,edgecolor =patch_colors[cindx],facecolor=patch_colors[cindx],alpha = patch_alpha))
        if patch_text: plt.text(0.302,2.24,'Protein',fontsize=11,alpha=1, color = 'k')
        cindx += 1
    #add on chemical reaction lines 
    #slopes are taken from Hatcher et al. (2003) Graphical method for analysis...
    #if 'hydrogenation' in plot_reactions: 
        # need to explicitly call the axes handle, and then add a rectangle specifying the (bottom left position), width, height
        #plt.plot((0.5,0.5),(2.0,0.5), "--",alpha=1,color='#d7191c')
        #plt.text(0.52,0.52,'H',fontsize=10,alpha=1,color='#d7191c')
    #if 'redox' in plot_reactions:
        #plt.plot((0.1,0.8),(1,1), "--",alpha=1,color='#fdae61')
        #plt.text(0.81,1.02,'R/O',fontsize=10,alpha=1,color='#fdae61')
    #if 'condensation' in plot_reactions:
        #plt.plot((0.2,0.8),(0.4,1.6), "--",alpha=1,color='#ffffbf')
        #plt.text(0.82,1.58,'Co',fontsize=10,alpha=1,color='#ffffbf')
    #if 'methylation' in plot_reactions: 
        #plt.plot((0.1,.6),(1.8,.8),"--",alpha = 1,color='#abd9e9')
        #plt.text(.61,.79,'Me',fontsize=10,alpha= 1,color='#abd9e9')
    #if 'carboxylation' in plot_reactions:
        #plt.plot((0.1,0.8),(2,2),"--",alpha=1,color='#2c7bb6')
        #plt.text(0.82,2.02,'Cbx',fontsize=10,alpha=1,color='#2c7bb6')
    plt.xlabel(f"Atomic ratio of {x_ratio[0]}/{x_ratio[1]}")
    plt.ylabel(f"Atomic ratio of {y_ratio[0]}/{y_ratio[1]}")
    fig = plt.gcf()
    ax = plt.gca()
    return fig, ax 

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