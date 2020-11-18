from ..formula import *
from .relative_intensity import relative_intensity
import numpy as np

def diversity_indices (formula_list,intensity_list,mz_list = [],indices = ['r','GS','SW','C','O','NOSC','DBE','rAI','HC','OC'],verbose = True):
    
    """ 
	Docstring for function pyKrev.diversity_indices
	====================
	This function takes a list of formula and list of corresponding peak intensities
	calculates a variety of diversity estimates (richness, abundance-based and functional.)
    
	Use
	----
	diversity_indices(Y,X)
    
	Returns a dictionary containing the values listed in 'indices'.   
    
	Parameters
	----------
	Y: A list of molecular formula strings. 
	X: A list or numpy array of peak intensities. 
	mz_values: A list of mz_values (optional). 
	indices: a list of strings specifying the specific diversity indices to calculate, can include:
        'r' : molecular richness (i.e. number of molecular formula)
        'GS': Gini-simpson abundance based alpha diversity (species eveness)
        'SW': Shannon Wiener abundance based alpha diversity
        'C' : C number functional diversity
        'O' : O number functional diversity
        'N' : N number functional diversity
        'HC': HC ratio functional diversity   
        'OC': OC ratio  functional diversity   
        'NC': NC ratio functional diversity 
        'NOSC': Nominal oxidation state of C functional diversity    
        'DBE' : Double bound equivalent functional diversity          
        'rAI' : reformulated aromaticity index functional diversity 
        'mz'  : mz values functional diversity
    verbose:  boolean, print the results to screen or not 
    
	Info
	----------
	For a greater indepth description of 
	diversity indices on molecular formula see: 'Mentges, Andrea, et al. "Functional molecular diversity of marine dissolved 
	organic matter is reduced during degradation." 
	Frontiers in Marine Science 4 (2017): 194.' 
    Functional diversity is calculated using Rao's quadratic entropy. See: "Rao, R. (1982). Diversity and dissimilarity coefficients: a unified approach." 
        
    """      
    formula_set = set(formula_list)
    if len(formula_set) != len(formula_list):
        print('Warning: duplicates detected in formula list. Remove to avoid inaccuracies.')
    if 'mz' in indices: 
        assert len(mz_list) == len(formula_list), 'you must provide an mz list if to calculate mz functional diversity'
    
    #calculate element counts, element ratios and other feature types for functional diversity metric
    count_list = element_counts(formula_list)
    ratio_list = element_ratios(formula_list, ratios = ['HC','OC','NC'])
    C_list = []
    O_list = []
    N_list = []
    HC_list = []
    OC_list = []
    NC_list = []
    
    for i,z in zip(count_list,ratio_list):
        C_list.append(i['C'])
        HC_list.append(z['HC'])
        O_list.append(i['O'])
        OC_list.append(z['OC'])
        N_list.append(i['N'])
        NC_list.append(z['NC'])

        
    C_list = np.array(C_list)
    O_list = np.array(O_list)
    HC_list = np.array(HC_list)
    OC_list = np.array(OC_list)
    N_list = np.array(N_list)
    NC_list = np.array(NC_list)
    
    AI_list = aromaticity_index(formula_list,index_type='rAI')
    DBE_list = double_bond_equivalent(formula_list)
    NOSC_s = nominal_oxidation_state(formula_list)
        
    diversity_indices = dict()
    
    #Firstly normalise the intensity_list
    rel_abundance = np.array(relative_intensity(intensity_list))
    
    #Molecular richness is just the count of molecular formula
    if 'r' in indices: 
        D_r = len([ i for i in rel_abundance if i!= 0])
        if verbose == True:
            print('Molecular richness:',D_r,'\n')
        diversity_indices['D_r'] = D_r
    
    #Abundance based diversity 
    if verbose == True:
        print('Abundanced based diversity:')
        
    #Gini-Simpson index decribes the distribution of compounds across molecular formulas
    #The Gini simpson index ranges from 0 to 1, where larger indiate higher diversity with a maximum value of 1 - 1/N for a uniform distribution. 
    
    if 'GS' in indices: 
        D_a_GS = 1 - sum(rel_abundance**2)
        if verbose == True:
            print('Gini-Simpson Index:', D_a_GS)
        diversity_indices['D_a_GS'] = D_a_GS
    
    #The Shannon Wiener Index accounts for the evenness of individual MFs in each sample.
    if 'SW' in indices: 
        D_a_SW = -sum(rel_abundance * np.log(rel_abundance,where=(rel_abundance!=0))) #where is useful because it stops us computing the natural log of 0 which is undefined
        if verbose == True:
            print('Shannon-Wiener Index:', D_a_SW,'\n')
        diversity_indices['D_a_SW'] = D_a_SW

    
    #Rao's entropy uses 'absolute difference' of a chemical characteristic which can be defined as |x - y| or abs(x - y)
    if verbose == True:
        print('Functional based diversity:')
    
    #Calculate functional diversity based on C number, H/C ratio and oxidation state of C 
    #Include modified aromaticity index 
    
    D_f_C = 0
    D_f_O = 0
    D_f_N = 0
    
    D_f_HC = 0 
    D_f_OC = 0
    D_f_NC = 0

    D_f_rAI = 0 
    D_f_DBE = 0 
    D_f_NOSC = 0 
    D_f_mz = 0 

    
    for i in range(0,len(formula_list)-1):
        if 'C' in indices: 
            D_f_C += sum((rel_abundance[i] * rel_abundance[i+1::] * abs(C_list[i] - C_list[i+1::])))
        if 'O' in indices: 
            D_f_O += sum((rel_abundance[i] * rel_abundance[i+1::] * abs(O_list[i] - O_list[i+1::])))
        if 'N' in indices: 
            D_f_N += sum(rel_abundance[i] * rel_abundance[i+1::] * abs(N_list[i] - N_list[i+1::]))
        if 'HC' in indices: 
            D_f_HC += sum(rel_abundance[i] * rel_abundance[i+1::] * abs(HC_list[i] - HC_list[i+1::]))
        if 'OC' in indices: 
            D_f_OC += sum(rel_abundance[i] * rel_abundance[i+1::] * abs(OC_list[i] - OC_list[i+1::]))
        if 'NC' in indices: 
            D_f_NC += sum(rel_abundance[i] * rel_abundance[i+1::] * abs(NC_list[i] - NC_list[i+1::]))
        if 'rAI' in indices:
            D_f_rAI += sum(rel_abundance[i] * rel_abundance[i+1::] * abs(AI_list[i] - AI_list[i+1::]))
        if 'DBE' in indices:
            D_f_DBE += sum(rel_abundance[i] * rel_abundance[i+1::] * abs(DBE_list[i] - DBE_list[i+1::]))
        if 'NOSC' in indices: 
            D_f_NOSC += sum(rel_abundance[i] * rel_abundance[i+1::] * abs(NOSC_s[i] - NOSC_s[i+1::]))
        if 'mz' in indices: 
            D_f_mz += sum(rel_abundance[i] * rel_abundance[i+1::] * abs(mz_list[i] - mz_list[i+1::]))
        
    
    if 'C' in indices: 
        if verbose == True:
            print('Raos Quadratic Index (C Number): ', D_f_C)
        diversity_indices['D_f_C'] = D_f_C
    
    if 'O' in indices: 
        if verbose == True:
            print('Raos Quadratic Index (O Number): ', D_f_O)
        diversity_indices['D_f_O'] = D_f_O
    
    if 'N' in indices: 
        if verbose == True:
            print('Raos Quadratic Index (N Number): ', D_f_N)
        diversity_indices['D_f_N'] = D_f_N
    
    if 'HC' in indices: 
        if verbose == True:
            print('Raos Quadratic Index (HC Ratio): ', D_f_HC)
        diversity_indices['D_f_HC'] = D_f_HC
    
    if 'OC' in indices:
        if verbose == True:
            print('Raos Quadratic Index (OC Ratio): ', D_f_OC)
        diversity_indices['D_f_OC'] = D_f_OC
    if 'NC' in indices: 
        if verbose == True:
            print('Raos Quadratic Index (NC Ratio): ', D_f_NC)
        diversity_indices['D_f_NC'] = D_f_NC
    
    if 'NOSC' in indices: 
        if verbose == True:
            print('Raos Quadratic Index (NOSC): ', D_f_NOSC)
        diversity_indices['D_f_NOSC'] = D_f_NOSC
    
    if 'rAI' in indices: 
        if verbose == True:
            print('Raos Quadratic Index (rAI): ', D_f_rAI)
        diversity_indices['D_f_rAI'] = D_f_rAI
    
    if 'DBE' in indices: 
        if verbose == True:
            print('Raos Quadratic Index (DBE): ', D_f_DBE)
        diversity_indices['D_f_DBE'] = D_f_DBE
    
    if 'mz' in indices: 
        if verbose == True:
            print('Raos Quadratic Index (mz): ', D_f_mz)
        diversity_indices['D_f_mz'] = D_f_mz
    
    return diversity_indices