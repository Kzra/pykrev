import unittest
import numpy as np
from pykrev import element_counts, element_ratios, double_bond_equivalent, aromaticity_index, nominal_oxidation_state, calculate_mass, kendrick_mass_defect, find_intersections, filter_spectral_interference, msTupleDict, msTuple, average_mstuple

class TestFORMULA(unittest.TestCase):

    def setUp(self):
        pass

    def test_element_counts(self):
        x = msTuple(['C13H14O5','C13H14N2O4S2','C36H45ClN6O12'],[],[])
        correct = [{'C':13,'H':14,'N':0,'O':5,'P':0,'S':0},{'C':13,'H':14,'N':2,'O':4,'P':0,'S':2},{'C':36,'H':45,'N':6,'O':12,'P':0,'S':0}]
        res = element_counts(x)
        self.assertEqual(res, correct)

    def test_element_ratios(self):
        x = msTuple(['C13H14O5','C13H14N2O4S2','C36H45ClN6O12'],[],[])
        correct = [{'HC':1.0769230769230769,'OC':0.38461538461538464,'NC':0.0},{'HC':1.0769230769230769,'OC':0.3076923076923077,'NC':0.15384615384615385},{'HC':1.25,'OC':0.3333333333333333,'NC':0.16666666666666666}]
        res = element_ratios(x, ratios = ['HC','OC','NC'])
        self.assertEqual(res, correct)    

    def test_double_bond_equivalent(self):
        x = msTuple(['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C32H34O8'],[],[])
        correct = np.array([7,8,17.5,16])
        res = double_bond_equivalent(x)
        self.assertIsNone(np.testing.assert_array_equal(res, correct))# returns None if the arrays are equal

    def test_ai(self):
        x = msTuple(['C32H34O8','C26H28O10','C12H6O8'],[],[])
        correct = np.array([0.33,0.19,0.5])
        res = aromaticity_index(x,index_type = 'AI')
        self.assertIsNone(np.testing.assert_array_equal(np.round(res,2), np.round(correct,2)))

    def test_rai(self):
        x = msTuple(['C11H11N1','C11H10O1','C27H14O6'],[],[])
        correct = np.array([0.545,0.545,0.56])
        res = aromaticity_index(x,index_type = 'rAI')
        self.assertIsNone(np.testing.assert_array_equal(np.round(res,2),np.round(correct,2)))

    def test_nosc(self):
        x = msTuple(['C3H7N1O2','C5H11N1O2','C5H11N1O2S1'],[],[])
        correct = np.array([0.00,-0.80,-0.40])
        res = nominal_oxidation_state(x)
        self.assertIsNone(np.testing.assert_array_equal(np.round(res,3),np.round(correct,3)))

    def test_mass_mono(self):
        x = msTuple(['C3H7N1O2','C5H11N1O2','C5H11N1O2P1S2'],[],[])
        correct = np.array([89.04768,117.0790,211.9969])
        res = calculate_mass(x,method = 'monoisotopic')
        self.assertIsNone(np.testing.assert_array_equal(np.round(res,3),np.round(correct,3)))

    def test_mass_mono_analyte(self):
        x = msTuple(['C46H43N1O35','C41H25N1O39','C36H62O40'],[],[])
        correct = np.array([1168.154286,1153.993093,1133.274460])
        res = calculate_mass(x,method = 'monoisotopic', protonated = True, ion_charge = -1)
        self.assertIsNone(np.testing.assert_array_equal(np.round(res,3),np.round(correct,3)))

    def test_mass_nominal(self):
        x = msTuple(['C3H7N1O2','C5H11N1O2','C5H11N1O2P1S2'],[],[])
        correct = np.array([89,117,212])
        res = calculate_mass(x,method = 'nominal')
        self.assertIsNone(np.testing.assert_array_equal(np.round(res,3),np.round(correct,3)))

    def test_mass_average(self):
        x = msTuple(['C3H7N1O2','C5H11N1O2','C5H11N1O2P1S2'],[],[])
        correct = np.array([89.09332,117.1466,212.252512])
        res =calculate_mass(x,method = 'average')
        self.assertIsNone(np.testing.assert_array_equal(np.round(res,3),np.round(correct,3)))

    def test_kendrick_mass(self):
        z = msTuple([],[],np.array([351.3269,449.4364,505.4989]))
        correct = np.array([350.9346,448.9346,504.9345])
        res, res2 =kendrick_mass_defect(z, base = ['CH2'])
        self.assertIsNone(np.testing.assert_array_equal(np.round(res,3),np.round(correct,3)))

    def test_kendrick_mass_defect(self):
        z = msTuple([],[],np.array([351.3269,449.4364,505.4989]))
        correct = np.array([0.065,0.065,0.066])
        res, res2 =kendrick_mass_defect(z,base = ['CH2'], rounding = 'even')
        self.assertIsNone(np.testing.assert_array_equal(np.round(res2,3),np.round(correct,3)))

    def test_kendrick_mass_defect_integer(self):
        z = msTuple([],[],np.array([351.3269,449.4364,505.4989]))
        correct = np.array([0.065,0.065,0.066])
        res, res2 =kendrick_mass_defect(z,base = ['CH2'], rounding = 'rint')
        self.assertIsNone(np.testing.assert_array_equal(np.round(res2,3),np.round(correct,3)))

    def test_kendrick_mass_defect_floor(self):
        z = msTuple([],[],np.array([351.3269,449.4364,505.4989]))
        correct = np.array([-0.935,-0.935,-0.934])
        res, res2 =kendrick_mass_defect(z,base = ['CH2'], rounding = 'floor')
        self.assertIsNone(np.testing.assert_array_equal(np.round(res2,3),np.round(correct,3)))

    def test_kendrick_mass_defect_ceil(self):
        z = msTuple([],[],np.array([351.3269,449.4364,505.4989]))
        correct = np.array([0.065,0.065,0.066])
        res, res2 =kendrick_mass_defect(z,base = ['CH2'], rounding = 'ceil')
        self.assertIsNone(np.testing.assert_array_equal(np.round(res2,3),np.round(correct,3)))

    def test_find_intersections(self):
        x = msTuple(['A','B','C','D'],np.array([1,2,3,4]),np.array([1,2,3,4]))
        x2 = msTuple(['A','B','D','E','F'],np.array([1,2,3,4,5]),np.array([1,2,3,4,5]))
        x3 = msTuple(['A','D','E','F'],np.array([1,2,3,4]),np.array([1,2,3,4]))
        R = msTupleDict()
        R['x'] = x
        R['x2'] = x2
        R['x3'] = x3
        res = find_intersections(R)
        self.assertEqual(res[('x','x2')],{'B'})
    
    def test_filter_si(self):
        x3 = np.array([98.4096,98.8121,136.2304])
        x2 = ["","",""]
        x = np.array([14982198,1195016,1039854])
        xt = msTuple(x2,x,x3)
        res1, res2, res3 = filter_spectral_interference(xt)
        self.assertIsNone(np.testing.assert_array_equal(res3,np.array([98.8121,136.2304])))
    
    def test_average_mstuple(self):
        x = msTuple(['C4H5O6','C5H6O7'], np.array([4,5]), np.array([120,5]))
        y = msTuple(['C4H5O6','C5H6O7'], np.array([6,6]), np.array([110,110]))
        z = msTuple(['C4H5O6','C5H6O7'], np.array([20,7]), np.array([100,110]))
        testDict = msTupleDict()
        testDict['x'] = x 
        testDict['y'] = y
        testDict['z'] = z
        #while average_mstuple occassionally returns the formula in a different order, the ordering of all the variables stays consistent
        averageTuple, stdStats = average_mstuple(testDict, intensityMethod = 'mean', mzMethod = 'mean', stdDev = True)
       
if __name__ == '__main__':
    unittest.main()