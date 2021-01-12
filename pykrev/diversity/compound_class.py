from ..formula import element_ratios 
from ..formula import aromaticity_index
from ..formula import element_counts
import pandas as pd
import os

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
        'FORM' - Compound classification based on the formularity algorithm. 
        'KEGG_BioMol' - Match molecular formula to those listed in 'compounds with biological roles' KEGG BRITE Pathways (databases located in ./compound_data).
        'KEGG_Phyto' - Match molecular formula to those listed in the 'phytochemicals' KEGG BRITE Heirarchy (databases located in ./compound_data). 
        'KEGG_Lipid' - Match molecular formula to those listed in the 'lipids' KEGG BRITE Heirarchy (databases located in ./compound_data). 
        'KEGG_Pesticide' - Match molecular formula to those listed in the 'pesticides' KEGG BRITE Heirarchy (databases located in ./compound_data). 
        'KEGG_Toxin' - Match molecular formula to those listed in the 'toxin' KEGG BRITE Heirarchy (databases located in ./compound_data). 
        'KEGG_All' - Match molecular formula to all of the possible categories.
        
    Notes: Please refer to the KEGG website for information on BRITE heirarchies. https://www.genome.jp/kegg/brite.html
        
    """
    count_list = element_counts(formula_list)
    cclassCounts = dict()
    compound_class = []
    
    #Retrieve the path to the diversity directory that this script is in, needed to retrieve compound class databases
    functionPath = os.path.realpath(__file__)
    functionDirectory = os.path.dirname(functionPath)
    
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
        aindex = aromaticity_index(formula_list, index_type = 'rAI')
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
    
    if method == 'FORM':
        cRatios = element_ratios(formula_list)
        cclassCounts['Lipid-like'] = 0
        cclassCounts['Carbohydrate-like'] = 0
        cclassCounts['Unsaturated hydrocarbons'] = 0
        cclassCounts['Condensed aromatics'] = 0
        cclassCounts['Lignin-like'] = 0
        cclassCounts['Tannin-like'] = 0
        cclassCounts['Amino Sugar-like'] = 0
        cclassCounts['Protein-like'] = 0
        cclassCounts['Not assigned'] = 0
        for ratio in cRatios:
            if ratio['OC'] >= 0.01 and ratio['OC'] <= 0.3 and ratio['HC'] >= 1.5 and ratio['HC'] <= 2.2:
                compound_class.append('Lipid-like')
                cclassCounts['Lipid-like'] += 1 
            elif ratio['OC'] >= 0.7 and ratio['OC'] <= 1.1 and ratio['HC'] >= 1.5 and ratio['HC'] <= 2.3:
                compound_class.append('Carbohydrate-like')
                cclassCounts['Carbohydrate-like'] += 1 
            elif ratio['OC'] >= 0.01 and ratio['OC'] <= 0.1 and ratio['HC'] >= 0.8 and ratio['HC'] <= 1.5:
                compound_class.append('Unsaturated hydrocarbons')
                cclassCounts['Unsaturated hydrocarbons'] += 1
            elif ratio['OC'] >= 0.01 and ratio['OC'] <= 1 and ratio['HC'] >= 0.2 and ratio['HC'] <= 0.8:
                compound_class.append('Condensed aromatics')
                cclassCounts['Condensed aromatics'] += 1
            elif ratio['OC'] >= 0.1 and ratio['OC'] <= 0.7 and ratio['HC'] >= 0.8 and ratio['HC'] <= 1.6:
                compound_class.append('Lignin-like')
                cclassCounts['Lignin-like'] += 1
            elif ratio['OC'] >= 0.7 and ratio['OC'] <= 1.2 and ratio['HC'] >= 0.8 and ratio['HC'] <= 1.6:
                compound_class.append('Tannin-like')
                cclassCounts['Tannin-like'] += 1
            elif ratio['OC'] >= 0.6 and ratio['OC'] <= 0.7 and ratio['HC'] >= 1.5 and ratio['HC'] <= 2.2:
                compound_class.append('Amino Sugar-like')
                cclassCounts['Tannin-like'] += 1
            elif ratio['OC'] >= 0.3 and ratio['OC'] <= 0.6 and ratio['HC'] >= 1.5 and ratio['HC'] <= 2.3:
                compound_class.append('Protein-like')
                cclassCounts['Protein-like'] += 1
        
    if 'KEGG' in method:
        if method == 'KEGG_BioMol':
            BRITE = pd.read_csv(f"{functionDirectory}\compound_data\Brite_BioMol_DF.csv")
        elif method == 'KEGG_Lipid': 
            BRITE = pd.read_csv(f"{functionDirectory}\compound_data\Brite_Lipid_DF.csv")
        elif method == 'KEGG_Phyto':
            BRITE = pd.read_csv(f"{functionDirectory}\compound_data\Brite_Phyto_DF.csv")
        elif method == 'KEGG_Pesticide': 
            BRITE = pd.read_csv(f"{functionDirectory}\compound_data\Brite_Pesticide_DF.csv")
        elif method == 'KEGG_Toxin': 
            BRITE = pd.read_csv(f"{functionDirectory}\compound_data\Brite_Toxin_DF.csv")
        elif method == 'KEGG_All': 
            BRITE = pd.read_csv(f"{functionDirectory}\compound_data\Brite_All_DF.csv")
        else:
            print('Error: KEGG Database Method not recognised. Refer to docstring.')
        BriteFormula = BRITE['F'].to_list()
        BriteCatA = BRITE['A'].to_list()
        KEGGCats = set(BriteCatA)
        cclassCounts['Not Matched'] = 0
        for c in KEGGCats:
            cclassCounts[c] = 0
        for f in formula_list:
            try:
                idx = BriteFormula.index(f)
                compound_class.append(BriteCatA[idx])
                cclassCounts[BriteCatA[idx]] += 1
            except ValueError:
                cclassCounts['Not Matched'] += 1
                compound_class.append('Not Matched')

    return compound_class,cclassCounts