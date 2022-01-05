import unittest
import numpy as np
from pykrev import diversity_indices, ordination_matrix, bray_curtis_matrix, compound_class, normalise_intensity, page_rank, msTuple, msTupleDict

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
        y = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        x = (y,z,[])
        correct = {'D_r':10}
        res = diversity_indices(x, indices = ['r'])
        self.assertEqual(res, correct)

    def test_GS(self):
        y = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        x = (y,z,[])
        correct = {'D_a_GS':0.8593}
        res = diversity_indices(x, indices = ['GS'])
        self.assertEqual(np.around(res['D_a_GS'],3), np.around(correct['D_a_GS'],3))

    def test_SW(self):
        y = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        x = (y,z,[])
        correct = {'D_a_SW':2.09}
        res = diversity_indices(x, indices = ['SW'])
        self.assertEqual(np.around(res['D_a_SW'],3), np.around(correct['D_a_SW'],3))

    def test_functionalC(self):
        y = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        x = (y,z,[])
        res = diversity_indices(x, indices = ['N'])

    def test_functionalNC(self):
        y = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        x = (y,z,[])
        res = diversity_indices(x, indices = ['NC'])

    def test_functionalrAI(self):
        y = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        x = (y,z,[])
        res = diversity_indices(x, indices = ['rAI'])

    def test_functionalmz(self):
        y = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        mz = np.array([232,340,132,904,321,431,3424,200,3204,1000])
        x = (y,z,mz)
        res = diversity_indices(x, indices = ['mz'])

    def test_ordination_matrix(self):
        x = msTuple(['A','B','C','D'],np.array([1,2,3,4]),np.array([1,2,3,4]))
        x2 = msTuple(['A','B','D','E','F'],np.array([1,2,3,4,5]),np.array([1,2,3,4,5]))
        x3 = msTuple(['A','D','E','F'],np.array([1,2,3,4]),np.array([1,2,3,4]))
        R = msTupleDict()
        ores = ordination_matrix(R)

    def test_normalise_ordination(self):
        x = msTuple(['A','B','C','D'],np.array([1,2,3,4]),np.array([1,2,3,4]))
        x2 = msTuple(['A','B','D','E','F'],np.array([1,2,3,4,5]),np.array([1,2,3,4,5]))
        x3 = msTuple(['A','D','E','F'],np.array([1,2,3,4]),np.array([1,2,3,4]))
        R = msTupleDict()
        ores = ordination_matrix(R)
        #n1res = normalise_intensity(ores)
        #n2res = normalise_intensity(ores, norm_subset = 'PPP', norm_method = 'binary')
        #n3res = normalise_intensity(ores, norm_subset = 'LOS', p_L = 3, norm_method = 'minmax')
        #n4res = normalise_intensity(ores, norm_subset = 'PPP', p_P = 0.73, norm_method = 'zscore')
        #n5res = normalise_intensity(ores, norm_subset = 'PPP', p_P = 0.02, norm_method = 'mean')
        #n6res = normalise_intensity(ores, norm_subset = 'LOS', p_P = 0.02, norm_method = 'mean', p_L = 1000)
        #n7res = normalise_intensity(ores, norm_subset = 'LOS', p_P = 0.02, norm_method = 'mean', p_L = 1000, log = True)
        #n8res = normalise_intensity(ores, norm_subset = 'ALL', p_P = 0.02, norm_method = 'none', p_L = 1000, log = True)
        
    def test_bray_curtis_matrix(self):
        x = msTuple(['A','B','C','D'],np.array([1,2,3,4]),np.array([1,2,3,4]))
        x2 = msTuple(['A','B','D','E','F'],np.array([1,2,3,4,5]),np.array([1,2,3,4,5]))
        x3 = msTuple(['A','D','E','F'],np.array([1,2,3,4]),np.array([1,2,3,4]))
        R = msTupleDict()
        ores = ordination_matrix(R)
        bres = bray_curtis_matrix(np.array(ores))

    def test_compound_class_MSCC(self):
        y = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        x = (y,[],z)
        res = compound_class(x, method = 'MSCC')

    def test_compound_class_KELL(self):
        y = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        x = (y,[],z)
        res = compound_class(x, method = 'KELL')

    def test_compound_class_FORM(self):
        y = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        x = (y,[],z)
        res = compound_class(x, method = 'FORM')

    def test_compound_class_KEGG(self):
        y = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S']
        z = np.array([1000,2432,3000,4201,2000,5990,1000,6520,8000,9001])
        x = msTuple(y,[],z)
        res = compound_class(x, method = 'KEGG_All')
     
    def test_page_rank(self):
        x = (['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C9H11NO2', 'C9H11NO3', 'C11H12N2O2', 'C5H7NO3', 'C5H9NO3', 'C6H12N2O4S2','C6H11NO3S'],[],[])
        correct = np.array([ 2.17651119,  2.17651119,  2.17651119, 21.73523322, 21.73523322, 2.17651119, 21.73523322, 21.73523322,  2.17651119,  2.17651119])
        res = page_rank(x)
        self.assertIsNone(np.testing.assert_array_equal(np.round(res,3),np.round(correct,3)))

if __name__ == '__main__':
    unittest.main()
