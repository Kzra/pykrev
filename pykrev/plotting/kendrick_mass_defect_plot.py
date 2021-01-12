from ..formula import kendrick_mass_defect
from matplotlib import pyplot as plt
def kendrick_mass_defect_plot(formula_list,mz_list, base = 'CH2', **kwargs):
    """ 
	Docstring for function PyKrev.kendrick_mass_defect_plot.py
	====================
	This function takes a list of molecular formula strings and a list or numpy array of mz values and creates a kendrick mass defect plot.
	Use
	----
	kendrick_mass_defect_plot(Y,X)
    
	Returns the figure and axes handles and a tuple containing two numpy arrays, the first contains the kendrick mass and the second the kendrick mass defect. 
    
	Parameters
	----------
	Y: A list of molecular formula. 
    X: A list or numpy array of mz_values of len(Y). 
	Base: Atom group used to define the Kendrick mass.
    **kwargs: key word arguments for pyplot.scatter(). 
    """
    
    kendrickMass, kendrickMassDefect = kendrick_mass_defect(formula_list,mz_list,base)
    plt.scatter(kendrickMass,kendrickMassDefect, **kwargs)
    plt.xlabel('Kendrick Mass')
    plt.ylabel('Kendrick Mass Defect')
    fig = plt.gcf()
    ax = plt.gca()
    return fig, ax, (kendrickMass, kendrickMassDefect)
    