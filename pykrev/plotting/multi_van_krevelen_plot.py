from matplotlib import pyplot as plt
from ..formula.element_ratios import element_ratios 
from ..formula.element_counts import element_counts
from matplotlib.patches import Rectangle
def multi_van_krevelen_plot(msTupleDict,
                            colours=[],
                            symbols=[], 
                            alphas = [],
                            edge_colours = [],
                            x_ratio = 'OC',
                            y_ratio = 'HC',
                            patch_classes = [],
                            patch_alpha = 0.3,
                            patch_colors = ['#762a83','#9970ab','#c2a5cf','#e7d4e8','#d9f0d3','#a6dba0','#5aae61','#1b7837'],
                            patch_text = True,
                            **kwargs):
    """ 
	Docstring for function pykrev.multi_van_krevelen_plot
	====================
	This function takes multiple msTuples and plots a van Krevelen diagram. 
    
	Use
	----
	multi_van_krevelen_plot(Y)
    
	Returns tuple containing figure and axes handles 
    
	Parameters
	----------
	Y: an msTupleDict
        
	colours: A list of strings providing marker colours corresponding to each msTuple in *Y. 
	See: https://matplotlib.org/2.1.1/api/_as_gen/matplotlib.pyplot.plot.html
    
	symbol: A list of scatter plot symbols to use corresponding to each msTuple in *Y. 
	See: https://matplotlib.org/3.2.1/api/markers_api.html
    
	edge_colours: a list of strings providing marker edge colours corresponding to each msTuple in *Y.
        
	alphas: list of floats providing transparency of the marker points corresponding to each msTuple in *Y.
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
    
    **kwargs: other key word arguments to pass to plt.scatter()
    """ 
    #Tests
    assert 'alpha' not in kwargs, 'provide a list of alpha values, alphas = ...'
    assert 'c' not in kwargs, 'provide a list of colour values, colours = ...'
    assert 'color' not in kwargs, 'provide a list of color values, colours = ...'
    assert 'marker' not in kwargs, 'provide a list of marker values, symbols = ...'
    assert 'edgecolors' not in kwargs, 'provide a list of edgecolors, edge_colours = ...'
    assert 'label' not in kwargs, 'provide a list of labels, group_labels = ...'
    #Setup
    ##apply colour blind safe colors taken from https://colorbrewer2.org/ 
    cols = ['#7fc97f','#beaed4','#fdc086','#74add1','#fdae61','#abd9e9','#fee090','#e0f3f8','#ffffbf']
    sybls = ['o','X','s','D','+','1','^','v','8']
    group_labels = list(msTupleDict.keys())
    if not colours:
        colours = [cols[i] for i in range(0,len(msTupleDict))]
    if not symbols: 
        symbols = [sybls[i] for i in range(0,len(msTupleDict))]
    if not alphas:
        alphas = [0.5] * len(msTupleDict)
    if not edge_colours:
        edge_colours = ['None'] * len(msTupleDict)
    assert len(msTupleDict) == len(group_labels) == len(colours) == len(symbols) == len(alphas) == len(edge_colours), 'Input variables must all be the same length'
    i = 0 
    for msTuple in msTupleDict.values(): 
        x_axis = []
        y_axis = []
        ratio_list = element_ratios(msTuple, ratios = [x_ratio,y_ratio])
        for ratios in ratio_list:
                x_axis.append(ratios[x_ratio])
                y_axis.append(ratios[y_ratio])
        plt.scatter(x_axis, y_axis, alpha=alphas[i], edgecolors=edge_colours[i],c=colours[i], marker = symbols[i], label = group_labels[i], **kwargs)
        i += 1 
    #apply grid lines 
    plt.grid(True) 
    #add on chemical class patches
    #boundaries taken from formularity software
    assert len(patch_colors) >= len(patch_classes), "Provide at least as many colors as classes"
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
    #label axis
    plt.xlabel(f"Atomic ratio of {x_ratio[0]}/{x_ratio[1]}")
    plt.ylabel(f"Atomic ratio of {y_ratio[0]}/{y_ratio[1]}")
    fig = plt.gcf()
    ax = plt.gca()
    return fig,ax