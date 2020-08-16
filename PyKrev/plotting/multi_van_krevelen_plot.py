from matplotlib import pyplot as plt

def multi_van_krevelen_plot(*ratio_lists,
                            group_labels = [],
                            colours=[],
                            symbols=[], 
                            alphas = [],
                            edge_colours = [],
                            **kwargs):
    
     
    """ 
	Docstring for function PyKrev.multi_van_krevelen_plot
	====================
	This function takes multiple lists of H/C-O/C ratios and plots a van Krevelen diagram. 
    
	Use
	----
	multi_van_krevelen_plot(*Y)
    
	Returns None 
    
	Parameters
	----------
	*Y: multiple ratio_lists (must contain H/C and O/C). See PyKrev.element_ratios.
    
	group_labels: A list of strings corresponding to each ratio_list in Y to be displayed in the legend. 
    
	colours: A list of strings providing marker colours corresponding to each ratio_list in Y. 
	See: https://matplotlib.org/2.1.1/api/_as_gen/matplotlib.pyplot.plot.html
    
	symbol: A list of scatter plot symbols to use corresponding to each ratio_list in Y. 
	See: https://matplotlib.org/3.2.1/api/markers_api.html
    
	edge_colours: a list of strings providing marker edge colours corresponding to each ratio_list in Y.
        
	alphas: list of floats providing transparency of the marker points corresponding to each ratio_list in Y. 
      
    """ 
    if len(ratio_lists) == 1: ratio_lists = ratio_lists[0] #this enables the user to pass a nested set of lists
    assert 'alpha' not in kwargs, 'provide a list of alpha values, alphas = ...'
    assert 'c' not in kwargs, 'provide a list of colour values, colours = ...'
    assert 'color' not in kwargs, 'provide a list of color values, colours = ...'
    assert 'marker' not in kwargs, 'provide a list of marker values, symbols = ...'
    assert 'edgecolors' not in kwargs, 'provide a list of edgecolors, edge_colours = ...'
    assert 'label' not in kwargs, 'provide a list of labels, group_labels = ...'
    
    cols = ['red','green','blue','purple','orange','magenta','cyan','yellow','black']
    sybls = ['o','X','s','D','+','1','^','v','8']
    
    
    if not group_labels: 
        group_labels = ['NA'] * len(ratio_lists)
    if not colours:
        colours = [cols[i] for i in range(0,len(ratio_lists))]
    if not symbols: 
        symbols = [sybls[i] for i in range(0,len(ratio_lists))]
    if not alphas:
        alphas = [0.5] * len(ratio_lists)
    if not edge_colours:
        edge_colours = ['None'] * len(ratio_lists)
    
    assert len(ratio_lists) == len(group_labels) == len(colours) == len(symbols) == len(alphas) == len(edge_colours), 'Input variables must all be the same length'

    i = 0 
    for ratio_list in ratio_lists: 

        x_axis = []
        y_axis = []

        for ratios in ratio_list:
                x_axis.append(ratios['OC'])
                y_axis.append(ratios['HC'])

        plt.scatter(x_axis, y_axis, alpha=alphas[i], edgecolors=edge_colours[i],c=colours[i], marker = symbols[i], label = group_labels[i], **kwargs)
        i += 1 
        
    #apply grid lines 
    plt.grid(True) 
    
    #label axis
    plt.xlabel('Atomic ratio of O/C')
    plt.ylabel('Atomic ratio of H/C')
    
    return 