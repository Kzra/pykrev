
from .multi_van_krevelen_plot import multi_van_krevelen_plot
from ..formula import missing_formula
from ..formula import element_counts
from ..formula import element_ratios 

def missing_plot(*groups,**kwargs):
    
    """ 
    
	Docstring for function PyKrev.missing_plot
	====================
	This function takes multiple lists of formula and plots a multi_van_Krevelen with the missing formula
	in each list represented by a different symbol. 
    
	Use
	----
	missing_plot(*Y)
    
	Returns a dictionary containing missing groups.   
    
	Parameters
	----------
	*Y: lists of atomic formula   
    
	**kwargs: key word arguments for multi_van_krevelen. See: pk.multi_van_krevelen_plot()
    
    """
    
    group_labels = []
    if 'group_labels' in kwargs: #if the user has supplied group_labels in kwargs, save them to a variable so unique_groups is labelled correctly 
        group_labels = kwargs['group_labels']
    
    
    missing_ratios = []
    missing_groups = missing_formula(*groups,group_labels = group_labels)
    for group in missing_groups:
        mc = element_counts(missing_groups[group])
        missing_ratios.append(element_ratios(mc))
   
    multi_van_krevelen_plot(missing_ratios,**kwargs) 
    
    return missing_groups
    