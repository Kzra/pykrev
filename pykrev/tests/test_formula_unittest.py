import unittest
import numpy as np
from pykrev import element_counts, element_ratios, double_bond_equivalent, aromaticity_index, nominal_oxidation_state, calculate_mass, kendrick_mass_defect, find_intersections, unique_formula, missing_formula, filter_spectral_interference

class TestFORMULA(unittest.TestCase):

    def setUp(self):
        pass

    def test_element_counts(self):
        x = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12']
        correct = [{'C':13,'H':14,'N':0,'O':5,'P':0,'S':0},{'C':13,'H':14,'N':2,'O':4,'P':0,'S':2},{'C':36,'H':45,'N':6,'O':12,'P':0,'S':0}]
        res = element_counts(x)
        self.assertEqual(res, correct)

    def test_element_ratios(self):
        x = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12']
        correct = [{'HC':1.0769230769230769,'OC':0.38461538461538464,'NC':0.0},{'HC':1.0769230769230769,'OC':0.3076923076923077,'NC':0.15384615384615385},{'HC':1.25,'OC':0.3333333333333333,'NC':0.16666666666666666}]
        res = element_ratios(x, ratios = ['HC','OC','NC'])
        self.assertEqual(res, correct)    

    def test_double_bond_equivalent(self):
        x = [{'C':13,'H':14,'N':0,'O':5,'P':0,'S':0},{'C':13,'H':14,'N':2,'O':4,'P':0,'S':2},{'C':36,'H':45,'N':6,'O':12,'P':0,'S':0},{'C':32,'H':34,'N':0,'O':8,'P':0,'S':0}]
        x = ['C13H14O5','C13H14N2O4S2','C36H45ClN6O12','C32H34O8']
        correct = np.array([7,8,17.5,16])
        res = double_bond_equivalent(x)
        self.assertIsNone(np.testing.assert_array_equal(res, correct))# returns None if the arrays are equal

    def test_ai(self):
        x = [{'C':32,'H':34,'N':0,'O':8,'P':0,'S':0},{'C':26,'H':28,'N':0,'O':10,'P':0,'S':0},{'C':12,'H':6,'N':0,'O':8,'P':0,'S':0}]
        x = ['C32H34O8','C26H28O10','C12H6O8']
        correct = np.array([0.33,0.19,0.5])
        res = aromaticity_index(x,index_type = 'AI')
        self.assertIsNone(np.testing.assert_array_equal(np.round(res,2), np.round(correct,2)))

    def test_rai(self):
        x = [{'C':11,'H':11,'N':1,'O':0,'P':0,'S':0},{'C':11,'H':10,'N':0,'O':1,'P':0,'S':0},{'C':27,'H':14,'N':0,'O':6,'P':0,'S':0}]
        x = ['C11H11N1','C11H10O1','C27H14O6']
        correct = np.array([0.545,0.545,0.56])
        res = aromaticity_index(x,index_type = 'rAI')
        self.assertIsNone(np.testing.assert_array_equal(np.round(res,2),np.round(correct,2)))

    def test_nosc(self):
        x = [{'C':3,'H':7,'N':1,'O':2,'P':0,'S':0},{'C':5,'H':11,'N':1,'O':2,'P':0,'S':0},{'C':5,'H':11,'N':1,'O':2,'P':0,'S':1}]
        x = ['C3H7N1O2','C5H11N1O2','C5H11N1O2S1']
        correct = np.array([0.00,-0.80,-0.40])
        res = nominal_oxidation_state(x)
        self.assertIsNone(np.testing.assert_array_equal(np.round(res,3),np.round(correct,3)))

    def test_mass_mono(self):
        x = ['C3H7N1O2','C5H11N1O2','C5H11N1O2P1S2']
        correct = np.array([89.04768,117.0790,211.9969])
        res = calculate_mass(x,method = 'monoisotopic')
        self.assertIsNone(np.testing.assert_array_equal(np.round(res,3),np.round(correct,3)))

    def test_mass_nominal(self):
        x = ['C3H7N1O2','C5H11N1O2','C5H11N1O2P1S2']
        correct = np.array([89,117,212])
        res = calculate_mass(x,method = 'nominal')
        self.assertIsNone(np.testing.assert_array_equal(np.round(res,3),np.round(correct,3)))

    def test_mass_average(self):
        x = ['C3H7N1O2','C5H11N1O2','C5H11N1O2P1S2']
        correct = np.array([89.09332,117.1466,212.252512])
        res =calculate_mass(x,method = 'average')
        self.assertIsNone(np.testing.assert_array_equal(np.round(res,3),np.round(correct,3)))

    def test_kendrick_mass(self):
        x = [ 'C23H43O2','C30H57O2','C34H65O2']
        z = np.array([351.3269,449.4364,505.4989])
        correct = np.array([350.9346,448.9346,504.9345])
        res, res2 =kendrick_mass_defect(x,mz_list = z,base = ['CH2'])
        self.assertIsNone(np.testing.assert_array_equal(np.round(res,3),np.round(correct,3)))

    def test_kendrick_mass_defect(self):
        x = [ 'C23H43O2','C30H57O2','C34H65O2']
        z = np.array([351.3269,449.4364,505.4989])
        correct = np.array([0.065,0.065,0.066])
        res, res2 =kendrick_mass_defect(x,mz_list = z,base = ['CH2'], rounding = 'even')
        self.assertIsNone(np.testing.assert_array_equal(np.round(res2,3),np.round(correct,3)))

    def test_kendrick_mass_defect_integer(self):
        x = [ 'C23H43O2','C30H57O2','C34H65O2']
        z = np.array([351.3269,449.4364,505.4989])
        correct = np.array([0.065,0.065,0.066])
        res, res2 =kendrick_mass_defect(x,mz_list = z,base = ['CH2'], rounding = 'integer')
        self.assertIsNone(np.testing.assert_array_equal(np.round(res2,3),np.round(correct,3)))

    def test_find_intersections(self):
        x = ['A','B','C','D']
        x2 = ['A','B','D','E','F']
        x3 = ['A','D','E','F']
        res = find_intersections(formula_lists =[x,x2,x3],group_labels = ['x','x2','x3'])
        self.assertEqual(res[('x',)],{'C'})
        self.assertEqual(res[('x','x2')],{'B'})

    def test_unique_formula(self):
        x = ['A','B','C','D']
        x2 = ['A','B','D','E','F']
        x3 = ['A','D','E','F']
        res = unique_formula(x,x2,x3, group_labels = ['x','x2','x3'])
        self.assertEqual(res['x'],['C'])
        self.assertEqual(res['x2'],[])

    def test_missing_formula(self):
        x = ['A','B','C','D']
        x2 = ['A','B','D','E','F']
        x3 = ['A','D','E','F']
        res = missing_formula(x,x2,x3, group_labels = ['x','x2','x3'])
        self.assertEqual(sorted(res['x']),sorted(['F','E']))
        self.assertEqual(res['x2'],[])
    
    def test_filter_si(self):
        x = [98.4096,98.8121,136.2304]
        x2 = ["","",""]
        x3 = [14982198,1195016,1039854]
        res1, res2, res3 = filter_spectral_interference(x,x2,x3)
        self.assertIsNone(np.testing.assert_array_equal(res1,np.array([98.8121,136.2304])))

        
if __name__ == '__main__':
    unittest.main()