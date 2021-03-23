import unittest
import numpy as np
from pykrev import diversity_indices, ordination_matrix, bray_curtis_matrix, compound_class, normalise_intensity

class TestDIVERSITY(unittest.TestCase):

    def setUp(self):
        pass

    def test_sum_relative_intensity(self):
        z = np.array([100,200,300])
        correct = np.array([100/600,200/600,300/600])
        res = normalise_intensity(z)
        self.assertIsNone(np.testing.assert_array_equal(np.round(res,3),np.round(correct,3)))

    def test_max_relative_intensity(self):
        z = np.array([100,200,300])
        correct = np.array([100/300,200/300,300/300])
        res = normalise_intensity(z, norm_method = 'max')
        self.assertIsNone(np.testing.assert_array_equal(np.round(res,3),np.round(correct,3)))

    def test_unit_relative_intensity(self):
        z = np.array([100,200,300])
        res = normalise_intensity(z, norm_method = 'unit_vector')
        self.assertEqual(np.round(sum(res**2),3),1.00)

    def test_zscore_relative_intensity(self):
        z = np.array([100,200,300,400,21,321,342,543])
        res = normalise_intensity(z, norm_method = 'zscore')
        self.assertEqual(np.round(np.mean(res),3),0.00)
        self.assertEqual(np.round(np.std(res),3),1.00)

    def test_minmax_relative_intensity(self):
        z = np.array([100,200,300,400,21,321,342,543])
        res = normalise_intensity(z, norm_method = 'minmax')
        self.assertEqual(min(res),0)
        self.assertEqual(max(res),1)

    def test_mean_relative_intensity(self):
        z = np.array([100,200,300,400,21,321,342,543])
        res = normalise_intensity(z, norm_method = 'mean')
        self.assertEqual(np.round(np.mean(res),3),0.00)
        
    def test_median_relative_intensity(self):
        z = np.array([100,200,300,400,21,321,342,543])
        res = normalise_intensity(z, norm_method = 'median')
        self.assertEqual(np.round(np.median(res),3),0.00)

    def test_binary_relative_intensity(self):
        z = np.array([100,0,300,400,21,321,0,543])
        res = normalise_intensity(z, norm_method = 'binary')
        self.assertEqual(sum(res),6)
        
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

    def test_normalise_ordination(self):
        x = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        x2 = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H31NO3', 'C11H12N1O2', 'C5H73O3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        z2 = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        ores = ordination_matrix(molecular_formulas = [x,x2],peak_intensities = [z,z2])
        n1res = normalise_intensity(ores)
        n2res = normalise_intensity(ores, norm_subset = 'PPP', norm_method = 'binary')
        n3res = normalise_intensity(ores, norm_subset = 'LOS', p_L = 3, norm_method = 'minmax')
        n4res = normalise_intensity(ores, norm_subset = 'PPP', p_P = 0.73, norm_method = 'zscore')
        n5res = normalise_intensity(ores, norm_subset = 'PPP', p_P = 0.02, norm_method = 'mean')
        n6res = normalise_intensity(ores, norm_subset = 'LOS', p_P = 0.02, norm_method = 'mean', p_L = 1000)
        n7res = normalise_intensity(ores, norm_subset = 'LOS', p_P = 0.02, norm_method = 'mean', p_L = 1000, log = True)
        n8res = normalise_intensity(ores, norm_subset = 'ALL', p_P = 0.02, norm_method = 'none', p_L = 1000, log = True)
        
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
        res = compound_class(x, method = 'KELL')

    def test_compound_class_FORM(self):
        x = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        res = compound_class(x, method = 'FORM')

    def test_compound_class_KEGG(self):
        x = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        res = compound_class(x, method = 'KEGG_All')

if __name__ == '__main__':
    unittest.main()
