
from .multi_van_krevelen_plot import multi_van_krevelen_plot
from ..formula import missing_formula
from ..formula import element_counts
from ..formula import element_ratios 

def missing_plot(*groups,group_labels=['NA','NA','NA','NA','NA'],colours=['blue','green','red','purple','black'], symbols=['o','^','X','D','s'], marker_size = 20, transparency = [0.5,0.5,0.5,0.5,0.5],
                        ref_compounds = ['lipid','condensed hydrocarbon', 'lignin', 'protein','cellulose'], 
                        ref_reactions = ['methylation','hydrogenation','hydration','redox'], edge_colours = ['none','none','none','none','none']):
    
    """ 
	Docstring for function PyKrev.missing_plot
	====================
	This function takes multiple lists of formula and plots a multi_van_Krevelen with the missing formula
	in each list represented by a different symbol. 
    
	Use
	----
	missing_plot(Y1,Y2,Y3,...,Yn)
    
	Returns figure, axes, and legend handles for a van krevelen plot and a dictionary containing missing groups.   
    
	Parameters
	----------
	Y1...Yn: lists of atomic formula   
    
	group_labels: A list of strings corresponding to each formula in Y to be displayed in the legend. 
    
	colours: A list of strings providing marker colours corresponding to each formula list in Y. 
	See: https://matplotlib.org/2.1.1/api/_as_gen/matplotlib.pyplot.plot.html
    
	symbol: A list of scatter plot symbols to use corresponding to each formula list in Y. 
	See: https://matplotlib.org/3.2.1/api/markers_api.html
    
	edge_colour: a list of strings providing marker edge colours corresponding to each formula list in Y.
    
	max_marker_size: list of floats providing marker size orresponding to each formula list in Y.
    
	trasparency: list of floats providing transparency of the marker points corresponding to each formula list in Y. 
      
	ref_compounds: a list of strings  
    
    """
    
    missing_ratios = []
    missing_groups = missing_formula(*groups)
    for group in missing_groups:
        mc = element_counts(missing_groups[group])
        missing_ratios.append(element_ratios(mc))
   
    fig, ax, legend = multi_van_krevelen_plot(ratio_lists = missing_ratios,group_labels = group_labels,colours=colours, symbols=symbols, marker_size = marker_size, transparency = transparency,
                        ref_compounds = ref_compounds)      
    
    return fig, ax, legend, missing_groups
    