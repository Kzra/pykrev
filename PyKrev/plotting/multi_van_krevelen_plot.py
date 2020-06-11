from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

def multi_van_krevelen_plot(ratio_lists=[],group_labels = ['NA','NA','NA','NA','NA'],colours=['blue','green','red','purple','black'],
                        symbols=['o','^','X','D','s'], marker_size = 20, transparency = [0.5,0.5,0.5,0.5,0.5],
                        ref_compounds = ['lipid','condensed hydrocarbon', 'lignin', 'protein','cellulose'], 
                        edge_colours = ['none','none','none','none','none'],subplot = None):
    
     
    """ 
	Docstring for function PyKrev.multi_van_krevelen_plot
	====================
	This function takes multiple lists of H/C-O/C ratios and plots a van Krevelen diagram. 
    
	Use
	----
	multi_van_krevelen_plot(Y1,Y2,Y3)
    
	Returns figure, axes, and legend handles for a van krevelen plot.   
    
	Parameters
	----------
	Y1..Yn: A list containing multiple nested ratio_lists (must contain H/C and O/C). See PyKrev.element_ratios.
    
	group_labels: A list of strings corresponding to each ratio_list in Y to be displayed in the legend. 
    
	colours: A list of strings providing marker colours corresponding to each ratio_list in Y. 
	See: https://matplotlib.org/2.1.1/api/_as_gen/matplotlib.pyplot.plot.html
    
	symbol: A list of scatter plot symbols to use corresponding to each ratio_list in Y. 
	See: https://matplotlib.org/3.2.1/api/markers_api.html
    
	edge_colour: a list of strings providing marker edge colours corresponding to each ratio_list in Y.
    
	max_marker_size: list of floats providing marker size orresponding to each ratio_list in Y.
    
	trasparency: list of floats providing transparency of the marker points corresponding to each ratio_list in Y. 
      
	ref_compounds: a list of strings  

             
    """ 
    if len(ratio_lists > 5):
        print('Warning: you must supply your own symbols, colours, group_labels and transparencies')
    
    legend_handles= []
    i = 0 
    
    #If this is being used by the subplot function it needs to inherit a figure and ax 
    if subplot == True:
        fig = plt.gcf()
        ax = plt.gca()
    #otherwise just create new fig and ax 
    else: 
        fig,ax = plt.subplots()
    
    for ratio_list in ratio_lists: 

        x_axis = []
        y_axis = []

        for ratios in ratio_list:
                x_axis.append(ratios['OC'])
                y_axis.append(ratios['HC'])

        scatter = ax.scatter(x_axis, y_axis, alpha=transparency[i], 
                             edgecolors=edge_colours[i], s=marker_size, c=colours[i], marker = symbols[i])
        legend_handles.append(Line2D([],[],marker=symbols[i],ls ="",alpha = transparency[i],mfc = colours[i],mec = edge_colours[i]))
        i += 1 
        
    #apply grid lines 
    ax.grid(True) 
    
    #set up legend
    legend = ax.legend(handles = legend_handles, labels = group_labels[0:i],loc="upper right")
 
    plt.xlabel('Atomic ratio of O/C')
    plt.ylabel('Atomic ratio of H/C')
    
    
    return fig, ax, legend 