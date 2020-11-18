import unittest
import numpy as np
from pykrev import van_krevelen_plot, element_ratios, element_counts, kendrick_mass_defect_plot, multi_van_krevelen_plot, van_krevelen_histogram, missing_plot, unique_plot, double_bond_equivalent 

class TestPLOTTING(unittest.TestCase):

    def setUp(self):
        pass

    def test_van_krevelen_plot(self):
       x = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
       x_dbe = double_bond_equivalent(x)
       van_krevelen_plot(x, y_ratio = 'NC', c = x_dbe)

    def test_multi_van_krevelen_plot(self):
       x = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
       x2 = ['C14H14O5','C12H14N2O4S2','C36H45ClN6O12','C9H14NO2', 'C9H11N2O3', 'C11H12N2O2', 'C55H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
       multi_van_krevelen_plot(x, x2, group_labels = ['Group_1','Group_2'])

    def test_missing_plot(self):
       x = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
       x2 = ['C14H14O5','C12H14N2O4S2','C36H45ClN6O12','C9H14NO2', 'C9H11N2O3', 'C11H12N2O2', 'C55H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
       missing_plot(x,x2,y_ratio = 'NC',group_labels = ['Group_1','Group_2'])

    def test_unique_plot(self):
       x = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
       x2 = ['C14H14O5','C12H14N2O4S2','C36H45ClN6O12','C9H14NO2', 'C9H11N2O3', 'C11H12N2O2', 'C55H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
       unique_plot(x,x2,y_ratio = 'SC',group_labels = ['Group_1','Group_2'])

    def test_kmd_plot(self):
       x = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
       z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
       kendrick_mass_defect_plot(x,z, base = ['CO'])

    def test_van_krevelen_histogram(self):
       x = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
       van_krevelen_histogram(x)

if __name__ == '__main__':
    unittest.main()