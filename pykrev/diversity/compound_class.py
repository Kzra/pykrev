from ..formula import element_ratios 
from ..formula import aromaticity_index
from ..formula import element_counts


def compound_class(formula_list, mass_list = [], method = 'MSCC'):
    """ 
	Docstring for function pyKrev.compound_class
	====================
	This function takes a list of formula and calculates compound class categories based on the method give in method. 
    
	Use
	----
	compound_class(Y)
    
	Returns two items. The first is a list of len(Y) in which each item is the compound class associated with the ith element of Y. 
    The second is a dictionary containing the compound class counts. 
    
	Parameters
	----------
	Y: A list of molecular formula strings 
    
	mass_list: A list of formula masses corresponding to the counts in element_counts. Required by MSCC method. 
    
	method: the method to use, one of:
        'MSCC' - Multi-dimensional stoiochiometric compound classification.  See: Albert Rivas-Ubach, Yina Liu, Thomas Stephen Bianchi, Nikola Tolic, Christer Jansson, and Ljiljana Paša-Tolic. (2018)
        Moving beyond the van Krevelen diagram: A new stoichiometric approach for compound classification in organisms. Analytical Chemistry, DOI: 10.1021/acs.analchem.8b00529
        'KELL' - Compound classification based on aromaticity index. See: Kellerman, A., Dittmar, T., Kothawala, D. et al. Chemodiversity of dissolved organic matter in lakes driven by climate and hydrology. Nat Commun 5, 3804 (2014).          https://doi.org/10.1038/ncomms4804
    """
    count_list = element_counts(formula_list)
    cclassCounts = dict()
    compound_class = []
    if method == 'MSCC':
        ## assert that the same number of elements in mass and count list
        assert len(mass_list) == len(count_list), 'to perform MSCC you must provide a mass list'
        ## initialise dictionary that contains counts 
        cclassCounts['Lipid'] = 0
        cclassCounts['Carbohydrate'] = 0
        cclassCounts['Amino-sugar'] = 0
        cclassCounts['Protein'] = 0
        cclassCounts['Phytochemical'] = 0
        cclassCounts['Nucleotide'] = 0
        cclassCounts['Not matched'] = 0
        cclassCounts['Double matched'] = 0
        ## initialise ratio list 
        cRatios = element_ratios(formula_list,ratios = ['OC','HC','NC','PC','NP'])
        for c,m,ratios in zip(count_list,mass_list,cRatios):
            ## a list that temporarily holds assigned categories
            cClass = []
            ## nucleotides cannot be double matched
            if ratios['OC'] >= 0.5 and ratios['OC'] < 1.7 and ratios['HC'] > 1 and ratios['HC'] < 1.8 and ratios['NC'] >= 0.2 and ratios['NC'] <= 0.5 and ratios['PC'] >= 0.1 and ratios['PC'] <= 0.35 and ratios['NP'] > 0.6 and ratios['NP'] <= 5 and c['N'] >= 2 and c['P'] >= 1 and c['S'] == 0 and m > 305 and m < 523:
                compound_class.append('Nucleotide')
                cclassCounts['Nucleotide'] += 1
            else:
                if ratios['OC'] <= 0.6 and ratios['HC'] >= 1.32 and ratios['NC'] <= 0.126 and ratios['PC'] < 0.35 and ratios['NP'] <= 5:
                    cClass.append('Lipid')
                if ratios['OC'] >= 0.8 and ratios['HC'] >= 1.65 and ratios['HC'] < 2.7  and c['N'] == 0:
                    cClass.append('Carbohydrate')
                if ratios['OC'] >= 0.61 and ratios['HC'] >= 1.45 and ratios['NC'] <= 0.2 and ratios['NC'] > 0.07 and ratios['PC'] < 0.3 and ratios['NP'] <= 2 and c['O'] >= 3 and c['N'] >= 1:
                    cClass.append('Amino-sugar')
                if ratios['OC'] <= 1.15 and ratios['HC'] < 1.32 and ratios['NC'] < 0.126 and ratios['PC'] <= 0.2 and ratios['NP'] <= 3:
                    cClass.append('Phytochemical')
                if ratios['OC'] > 0.12 and ratios['OC'] <= 0.6 and ratios['HC'] > 0.9 and ratios['HC'] < 2.5 and ratios['NC'] >= 0.126 and ratios['NC'] <= 0.7  and ratios['PC'] < 0.17 and c['N'] >= 1:
                    cClass.append('Protein')
                elif ratios['OC'] > 0.6 and ratios['OC'] <= 1 and ratios['HC'] > 1.2 and ratios['HC'] < 2.5 and ratios['NC'] > 0.2 and ratios['NC'] <= 0.7 and ratios['PC'] < 0.17 and c['N'] >= 1: 
                    cClass.append('Protein')
                if len(cClass) > 1: 
                    compound_class.append(f"Double matched: {cClass[0]} {cClass[1]}")
                    cclassCounts['Double matched'] += 1
                elif len(cClass) == 0:
                    compound_class.append("Not matched")
                    cclassCounts['Not matched'] += 1
                else:
                    compound_class.append(cClass[0])
                    cclassCounts[cClass[0]] += 1 
    if method == 'KELL':
        cRatios = element_ratios(formula_list)
        aindex = aromaticity_index(formula_list, index_type = 'AI')
        cclassCounts['Combustion-derived polycyclic aromatics'] = 0
        cclassCounts['Vascular plant-derived polyphenols'] = 0
        cclassCounts['Highly unsaturated and phenolic compounds'] = 0
        cclassCounts['Aliphatic compounds'] = 0
        for ai,ratio in zip(aindex,cRatios):
            if ai > 0.66: 
                compound_class.append('Combustion-derived polycyclic aromatics')
                cclassCounts['Combustion-derived polycyclic aromatics'] += 1 
            elif ai > 0.5 and ai <= 0.66:
                compound_class.append('Vascular plant-derived polyphenols')
                cclassCounts['Vascular plant-derived polyphenols'] += 1 
            elif ai <= 0.5 and ratio['HC'] < 1.5: 
                compound_class.append('Highly unsaturated and phenolic compounds')
                cclassCounts['Highly unsaturated and phenolic compounds'] += 1 
            elif ai <= 0.5 and ratio['HC'] >= 1.5:
                compound_class.append('Aliphatic compounds')
                cclassCounts['Aliphatic compounds'] += 1     
                
    return compound_class,cclassCounts