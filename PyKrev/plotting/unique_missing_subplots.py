from .multi_van_krevelen_plot import multi_van_krevelen_plot
from ..formula import unique_formula
from ..formula import missing_formula
from ..formula import element_counts
from ..formula import element_ratios 
import numpy as np
from matplotlib import pyplot as plt 

def unique_missing_subplots(*groups,titles = []):
    
    """ 
	Docstring for function PyKrev.unique_missing_subplots
	====================
	This function takes multiple lists of formula and plots multiple multi_van_krevelen diagrams
	with each plot containing the unique and missing formula associated with an input list.  
    
	Use
	----
	unique_missing_subplots(Y1,Y2,Y3,...,Yn)
    
	Returns figure, axes subplot handles and a unique_formula and missing_formula dictionary.  
    
	Parameters
	----------
	Y1...Yn: lists of atomic formula   
    
	titles: titles for each plot corresponding to each list in Y. 

	Example
	----------
    
	figs,axs,u,m = unique_missing_subplots(Dgroup2,Dgroup3,Dgroup4,Dgroup5,Dgroup6,Dgroup7,Dgroup8,Dgroup9,
	titles = ['A','B','C','D','E','F','G','H'])
   
    
    """
    
    n_groups = len(groups)
    
    if not titles: 
        titles = ['NA'] * n_groups
    else: 
        if len(titles) != n_groups:
            raise Exception('You need to supply as many titles as groups.')
            
    if n_groups < 2:
            raise Exception('You need to supply at least two groups.')
    
    i = 0
    n = 0 
    
    if n_groups > 2:
        columns = int(np.ceil(n_groups/2))
        tcol = columns - 1 
    else:
        columns = 2
        tcol = 2
        
    #figs and axs control the entire plot 
    figs, axs = plt.subplots( 2,columns,sharey=True) #the figures are shown in a 2 * n configuration
        
    missing_groups = missing_formula(*groups)
    unique_groups = unique_formula(*groups)
    
    for mgroup,ugroup,title in zip(missing_groups,unique_groups,titles):
        mc = element_counts(missing_groups[mgroup])
        uc = element_counts(unique_groups[ugroup])
        mr = element_ratios(mc)
        ur = element_ratios(uc)
        plt.sca(axs[n,i])
        #fig and ax controls each specific subplot
        fig,ax,legend = multi_van_krevelen_plot(ratio_lists = [mr,ur],group_labels = ['Missing','Unique'],symbols = ['X','o'],colours = ['r','g'],transparency = [0.8,0.5], subplot = True)
        #fig.set_size_inches(10,15)
        ax.text(0.05,3.75,title,fontsize = 12)
        
        if i == tcol:  #iterate through the subplots left to right
            n += 1 
            i = 0 
        else:
            i += 1 
    
    if n_groups == 2: # delete the unused subplots
       figs.delaxes(axs[1,0])
       figs.delaxes(axs[1,1])
    
    elif n_groups % 2 != 0:
        figs.delaxes(axs[1,tcol]) #if there is an odd number of figures to plot, delete the last subplot
        
    if n_groups <5: 
        figs.set_size_inches(20,15) #this ensures the plots scale nicely. 
    elif n_groups < 7:   
        figs.set_size_inches(25,15)
    else: 
        figs.set_size_inches(33,15)
        
        
    return figs,axs,unique_groups,missing_groups

