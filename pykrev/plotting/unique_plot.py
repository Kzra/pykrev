
from .multi_van_krevelen_plot import multi_van_krevelen_plot
from ..formula import unique_formula
from ..formula import element_counts
from ..formula import element_ratios 

def unique_plot(*groups, **kwargs):
    
    """ 
	Docstring for function PyKrev.unique_plot
	====================
	This function takes multiple lists of formula and plots a multi_van_Krevelen with the unique formula
	in each list represented by a different symbol. 
    
	Use
	----
	unique_plot(Y*)
    
	Returns figure and ax handles and the unique_groups produced by pk.unique_formula().   
    
	Parameters
	----------
	Y*: lists of atomic formula   
    
	**kwargs: key word arguments for pk.multi_van_krevelen_plot. See. pk.multi_van_krevelen_plot for more information.
    
    """ 
    #assert ('group_labels' in kwargs), 'you must provide group_labels'
    group_labels = []
    if 'group_labels' in kwargs: #if the user has supplied group_labels in kwargs, save them to a variable so unique_groups is labelled correctly 
        group_labels = kwargs['group_labels']
    
    unique_formulas = []
    unique_groups = unique_formula(*groups,group_labels = group_labels)
    
    for group in unique_groups:
       unique_formulas.append(unique_groups[group])
    
    fig,ax = multi_van_krevelen_plot(unique_formulas, **kwargs)
    
    return fig,ax,unique_groups
    