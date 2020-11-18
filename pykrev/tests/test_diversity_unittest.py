import unittest
import numpy as np
from pykrev import diversity_indices, relative_intensity, ordination_matrix, bray_curtis_matrix, compound_class 

class TestDIVERSITY(unittest.TestCase):

    def setUp(self):
        pass

    def test_relative_intensity(self):
        z = np.array([100,200,300])
        correct = np.array([100/600,200/600,300/600])
        res = relative_intensity(z)
        self.assertIsNone(np.testing.assert_array_equal(np.round(res,3),np.round(correct,3)))

    def test_richness(self):
        x = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        correct = {'D_r':10}
        res = diversity_indices(x,z, indices = ['r'])
        self.assertEqual(res, correct)

    def test_GS(self):
        x = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        correct = {'D_a_GS':0.8593}
        res = diversity_indices(x,z, indices = ['GS'])
        self.assertEqual(np.around(res['D_a_GS'],3), np.around(correct['D_a_GS'],3))

    def test_SW(self):
        x = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        correct = {'D_a_SW':2.09}
        res = diversity_indices(x,z, indices = ['SW'])
        self.assertEqual(np.around(res['D_a_SW'],3), np.around(correct['D_a_SW'],3))

    def test_functionalC(self):
        x = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        res = diversity_indices(x,z, indices = ['N'])

    def test_functionalNC(self):
        x = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        res = diversity_indices(x,z, indices = ['NC'])

    def test_functionalrAI(self):
        x = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        res = diversity_indices(x,z, indices = ['rAI'])

    def test_functionalmz(self):
        x = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        mz = np.array([232,340,132,904,321,431,3424,200,3204,1000])
        res = diversity_indices(x,z, mz_list = mz, indices = ['mz'])

    def test_ordination_matrix(self):
        x = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        x2 = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H31NO3', 'C11H12N1O2', 'C5H73O3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        z2 = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        ores = ordination_matrix(molecular_formulas = [x,x2],peak_intensities = [z,z2])

    def test_bray_curtis_matrix(self):
        x = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        x2 = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H31NO3', 'C11H12N1O2', 'C5H73O3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        z2 = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        ores = ordination_matrix(molecular_formulas = [x,x2],peak_intensities = [z,z2])
        bres = bray_curtis_matrix(np.array(ores))

    def test_compound_class_MSCC(self):
        x = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        res = compound_class(x,mass_list =z, method = 'MSCC')

    def test_compound_class_KELL(self):
        x = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        res = compound_class(x,mass_list=z, method = 'KELL')

if __name__ == '__main__':
    unittest.main()