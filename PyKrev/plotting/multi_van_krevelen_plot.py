from matplotlib import pyplot as plt
from ..formula import element_ratios 
from ..formula import element_counts

def multi_van_krevelen_plot(*formula_lists,
                            group_labels = [],
                            colours=[],
                            symbols=[], 
                            alphas = [],
                            edge_colours = [],
                            x_ratio = 'OC',
                            y_ratio = 'HC',
                            **kwargs):

    """ 
	Docstring for function PyKrev.multi_van_krevelen_plot
	====================
	This function takes multiple lists of molecula formula strings and plots a van Krevelen diagram. 
    
	Use
	----
	multi_van_krevelen_plot(*Y)
    
	Returns None 
    
	Parameters
	----------
	*Y: multiple lists of molecular formula strings.
    
	group_labels: A list of strings corresponding to each formula_list in Y to be displayed in the legend. 
    
	colours: A list of strings providing marker colours corresponding to each formula_list in Y. 
	See: https://matplotlib.org/2.1.1/api/_as_gen/matplotlib.pyplot.plot.html
    
	symbol: A list of scatter plot symbols to use corresponding to each formula_list in Y. 
	See: https://matplotlib.org/3.2.1/api/markers_api.html
    
	edge_colours: a list of strings providing marker edge colours corresponding to each formula_list in Y.
        
	alphas: list of floats providing transparency of the marker points corresponding to each formula_list in Y.
    x_ratio: element ratio to plot on x axis, given numerator denominator e.g. 'OC'
    y_ratio: element ratio to plot on y axis, given numerator denominator e.g. 'HC'
    **kwargs: other key word arguments to pass to plt.scatter()
      
    """ 
    if len(formula_lists) == 1: formula_lists = formula_lists[0] #this enables the user to pass a nested set of lists
    assert 'alpha' not in kwargs, 'provide a list of alpha values, alphas = ...'
    assert 'c' not in kwargs, 'provide a list of colour values, colours = ...'
    assert 'color' not in kwargs, 'provide a list of color values, colours = ...'
    assert 'marker' not in kwargs, 'provide a list of marker values, symbols = ...'
    assert 'edgecolors' not in kwargs, 'provide a list of edgecolors, edge_colours = ...'
    assert 'label' not in kwargs, 'provide a list of labels, group_labels = ...'
    
    cols = ['red','green','blue','purple','orange','magenta','cyan','yellow','black']
    sybls = ['o','X','s','D','+','1','^','v','8']
    
    
    if not group_labels: 
        group_labels = ['NA'] * len(formula_lists)
    if not colours:
        colours = [cols[i] for i in range(0,len(formula_lists))]
    if not symbols: 
        symbols = [sybls[i] for i in range(0,len(formula_lists))]
    if not alphas:
        alphas = [0.5] * len(formula_lists)
    if not edge_colours:
        edge_colours = ['None'] * len(formula_lists)
    
    assert len(formula_lists) == len(group_labels) == len(colours) == len(symbols) == len(alphas) == len(edge_colours), 'Input variables must all be the same length'

    i = 0 
    for formula_list in formula_lists: 

        x_axis = []
        y_axis = []

        ratio_list = element_ratios(formula_list, ratios = [x_ratio,y_ratio])

        for ratios in ratio_list:
                x_axis.append(ratios[x_ratio])
                y_axis.append(ratios[y_ratio])

        plt.scatter(x_axis, y_axis, alpha=alphas[i], edgecolors=edge_colours[i],c=colours[i], marker = symbols[i], label = group_labels[i], **kwargs)
        i += 1 
        
    #apply grid lines 
    plt.grid(True) 
    
    #label axis
    plt.xlabel(f"Atomic ratio of {x_ratio[0]}/{x_ratio[1]}")
    plt.ylabel(f"Atomic ratio of {y_ratio[0]}/{y_ratio[1]}")
    
    return 