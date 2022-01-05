import unittest
import numpy as np
from pykrev import van_krevelen_plot, element_ratios, element_counts, kendrick_mass_defect_plot, multi_van_krevelen_plot, van_krevelen_histogram, double_bond_equivalent, atomic_class_plot, compound_class_plot, mass_spectrum, mass_histogram, reaction_network, msTupleDict, msTuple

class TestPLOTTING(unittest.TestCase):

    def setUp(self):
        pass

    def test_van_krevelen_plot_no_patch(self):
       x = (['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S'],[],[])
       x_dbe = double_bond_equivalent(x)
       van_krevelen_plot(x, y_ratio = 'NC', c = x_dbe, patch_classes = [])

    def test_van_krevelen_plot_patch(self):
       x = (['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S'],[],[])
       x_dbe = double_bond_equivalent(x)
       fig,ax = van_krevelen_plot(x, y_ratio = 'NC', c = x_dbe, patch_classes = ['lipid-like','lignin-like'], patch_alpha = 0.2, patch_text = False, patch_colors = ['#ffffbf','#ffffbf'])

    def test_multi_van_krevelen_plot_no_patch(self):
      x = msTuple(['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S'],[],[])
      x2 = msTuple(['C14H14O5','C12H14N2O4S2','C36H45ClN6O12','C9H14NO2', 'C9H11N2O3', 'C11H12N2O2', 'C55H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S'],[],[])
      R = msTupleDict()
      R['x'] = x
      R['x2'] = x2
      multi_van_krevelen_plot(R, patch_classes = [])

    def test_multi_van_krevelen_plot_patch(self):
       x = msTuple(['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S'],[],[])
       x2 = msTuple(['C14H14O5','C12H14N2O4S2','C36H45ClN6O12','C9H14NO2', 'C9H11N2O3', 'C11H12N2O2', 'C55H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S'],[],[])
       R = msTupleDict()
       R['x'] = x
       R['x2'] = x2
       multi_van_krevelen_plot(R, patch_classes = ['tannin-like','lignin-like'], patch_alpha = 0.4, patch_text = True, patch_colors = ['#ffffbf','#ffffbf'])

    def test_kmd_plot(self):
       z = ([],[],np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001]))
       kendrick_mass_defect_plot(z, base = ['CO'], rounding = 'even')
       kendrick_mass_defect_plot(z, base = ['CO'], rounding = 'rint')
       kendrick_mass_defect_plot(z, base = ['CO'], rounding = 'ceil')
       kendrick_mass_defect_plot(z, base = ['CO'], rounding = 'floor')

    def test_van_krevelen_histogram(self):
       x = (['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S'],[],[])
       van_krevelen_histogram(x)

    def test_atomic_class_plot(self):
       x = (['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S'],[],[])
       atomic_class_plot(x, element = 'S')

    def test_compound_class_plot(self):
       x = (['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S'],[],[])
       compound_class_plot(x, method = 'KELL')
        
    def test_mass_histogram(self):
       x = (['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S'],[],[1,2,3,4,5,6,7,8,9,10])
       mass_histogram(x)
       mass_histogram(x, method = 'nominal')
    
    def test_mass_spectrum(self):
       y = np.array([3210,43,432,423,42,10,103,305,2054,1388])
       x = (['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S'],y,[])
       mass_spectrum(x)
       mass_spectrum(x, method = 'nominal')

    def test_reaction_network(self):
       y = np.array([3210,43,432,423,42,10,103,305,2054,1388])
       x = (['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S'],y,[])
       reaction_network(x)
       reaction_network(x, nodeAnnotations = {'Peak Intensity' : y})
   
if __name__ == '__main__':
    unittest.main()