from ..formula import *
from .relative_intensity import relative_intensity
import numpy as np

def diversity_indices (formula_list,intensity_list,indices = ['r','GS','SW','C','O','NOSC','DBE','rAI','HC','OC']):
    
    """ 
	Docstring for function pyKrev.diversity_indices
	====================
	This function takes a list of formula and list of corresponding peak intensities
	calculates a variety of diversity estimates (richness, abundance-based and functional.)
    
	Use
	----
	diversity_indices(Y,X)
    
	Returns a list of len(Y) in which each item is a dictionary containing the values listed in 'indices'.   
    
	Parameters
	----------
	Y: A list of molecular formula strings. 
	X: A list of peak intensities. 
	indices: a list of strings specifying the specific diversity indices to calculate:
	'r' : molecular richness
	'GS': Gini-simpson
	'SW': Shannon Wiener
	'C' : C number
	'O' : O number
	'HC': HC ratio   
	'OC': OC ratio    
	'NOSC': Nominal oxidation state of C    
	'DBE' : Double bound equivalence          
	'rAI' : modified aromaticity index 
    
	Info
	----------
	For a greater indepth description of 
	diversity indeces on molecular formula see: 'Mentges, Andrea, et al. "Functional molecular diversity of marine dissolved 
	organic matter is reduced during degradation." 
	Frontiers in Marine Science 4 (2017): 194.' 
        
    """      
    
    
    #calculate element counts, element ratios and other feature types for functional diversity metric
    count_list = element_counts(formula_list)
    ratio_list = element_ratios(count_list)
    C_list = []
    O_list = []
    HC_list = []
    OC_list = []
    
    for i,z in zip(count_list,ratio_list):
        C_list.append(i['C'])
        HC_list.append(z['HC'])
        O_list.append(i['O'])
        OC_list.append(z['OC'])
        
    C_list = np.array(C_list)
    O_list = np.array(O_list)
    HC_list = np.array(HC_list)
    OC_list = np.array(OC_list)
    AI_list = np.array(aromaticity_index(count_list,index_type='mod'))
    DBE_list = np.array(double_bond_equivalent(count_list))
    NOSC_s = np.array(nominal_oxidation_state(count_list))
        
    diversity_indices = dict()
    
    #Firstly normalise the intensity_list
    rel_abundance = np.array(relative_intensity(intensity_list))
    
    #Molecular richness is just the count of molecular formula
    if 'r' in indices: 
        D_r = len([ i for i in rel_abundance if i!= 0])
        print('Molecular richness:',D_r,'\n')
        diversity_indices['D_r'] = D_r
    
    #Abundance based diversity 
    print('Abundanced based diversity:')
        
    #Gini-Simpson index decribes the distribution of compounds across molecular formulas
    #The gene simpson index ranges from 0 to 1, where larger indiate higher diversity with a maximum value of 1 - 1/N for a uniform distribution. 
    
    if 'GS' in indices: 
        D_a_GS = 1 - sum(rel_abundance**2)
        print('Gini-Simpson Index:', D_a_GS)
        diversity_indices['D_a_GS'] = D_a_GS
    
    #The Shannon Wiener Index accounts for the evenness of individual MFs in each sample.
    if 'SW' in indices: 
        D_a_SW = -sum(rel_abundance * np.log(rel_abundance,where=(rel_abundance!=0))) #where is useful because it stops us computing the natural log of 0 which is undefined
        print('Shannon-Wiener Index:', D_a_SW,'\n')
        diversity_indices['D_a_SW'] = D_a_SW

    
    #Rao's entropy uses 'absolute difference' of a chemical characteristic which can be defined as |x - y| or abs(x - y)
    print('Functional based diversity:')
    
    #Calculate functional diversity based on C number, H/C ratio and oxidation state of C 
    #Include modified aromaticity index 
    
    D_f_C = 0
    D_f_O = 0
    
    D_f_NOSC = 0 
    
    D_f_HC = 0 
    D_f_OC = 0

    D_f_MOD_AI = 0 
    D_f_DBE = 0 

    
    for i in range(0,len(formula_list)-1):
        if 'C' in indices: 
            D_f_C += sum((rel_abundance[i] * rel_abundance[i+1::] * abs(C_list[i] - C_list[i+1::])))
        if 'O' in indices: 
            D_f_O += sum((rel_abundance[i] * rel_abundance[i+1::] * abs(O_list[i] - O_list[i+1::])))
        if 'NOSC' in indices: 
            D_f_NOSC += sum(rel_abundance[i] * rel_abundance[i+1::] * abs(NOSC_s[i] - NOSC_s[i+1::]))
        if 'HC' in indices: 
            D_f_HC += sum(rel_abundance[i] * rel_abundance[i+1::] * abs(HC_list[i] - HC_list[i+1::]))
        if 'OC' in indices: 
            D_f_OC += sum(rel_abundance[i] * rel_abundance[i+1::] * abs(OC_list[i] - OC_list[i+1::]))
        if 'rAI' in indices:
            D_f_MOD_AI += sum(rel_abundance[i] * rel_abundance[i+1::] * abs(AI_list[i] - AI_list[i+1::]))
        if 'DBE' in indices:
            D_f_DBE += sum(rel_abundance[i] * rel_abundance[i+1::] * abs(DBE_list[i] - DBE_list[i+1::]))
        
    
    if 'C' in indices: 
        print('Raos Quadratic Index (C Number): ', D_f_C)
        diversity_indices['D_f_C'] = D_f_C
    
    if 'O' in indices: 
        print('Raos Quadratic Index (O Number): ', D_f_O)
        diversity_indices['D_f_O'] = D_f_O
    
    if 'HC' in indices: 
        print('Raos Quadratic Index (HC Ratio): ', D_f_HC)
        diversity_indices['D_f_HC'] = D_f_HC
    
    if 'OC' in indices:
        print('Raos Quadratic Index (OC Ratio): ', D_f_OC)
        diversity_indices['D_f_OC'] = D_f_OC

    if 'NOSC' in indices:     
        print('Raos Quadratic Index (NOSC): ', D_f_NOSC)
        diversity_indices['D_f_NOSC'] = D_f_NOSC
    
    if 'rAI' in indices: 
        print('Raos Quadratic Index (MOD AI): ', D_f_MOD_AI)
        diversity_indices['D_f_MOD_AI'] = D_f_MOD_AI
    
    if 'DBE' in indices: 
        print('Raos Quadratic Index (DBE): ', D_f_DBE)
        diversity_indices['D_f_DBE'] = D_f_DBE
    
    return diversity_indices