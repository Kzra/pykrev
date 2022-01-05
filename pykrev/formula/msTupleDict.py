import pandas as pd
import numpy as np
from ..diversity.ordination_matrix import ordination_matrix
from .find_intersections import find_intersections

class msTupleDict(dict):
    """ 
    Docstring for class pykrev.msTupleDict
    ==========
    Core DataType for PyKrev batch processing. 
    A dictionary with sample names as keys and msTuple objects as values.

    Use
    ----------
    msTupleDict() 
    
    Returns an empty msTupleDict. 
    
    Methods 
    ----------
    msTupleDict.validate(): validate all of the msTuples stored in the dictionary.

    msTupleDict.summary(): summarise all of the msTuples stored in the dictionary.

    msTupleDict.intersections(exclusive = True): return a dictionary contanining all intersections between the formula in msTupleDict. See pk.find_intersections.

    msTupleDict.to_OrdinationMatrix(impute_value = 'nan'): write the contents of the msTupleDict to an ordination matrix. See pk.ordination_matrix. 

    msTupleDict.to_DataFrame(): write the contents of the msTupleDict to a pandas dataframe. Columns are 'assigned formula', 'mean mz' and 'std mz'

    """
    def validate(self):
        for v in self.values():
            v.validate()

    def summary(self):
        print(f'msTupleDict containing {len(self)} samples.')
        print()
        for k,v in zip(self.keys(),self.values()):
            print(f'{k} summary')
            print(f'{"*"*20}')
            v.summary()
            print()

    def intersections(self, exclusive = True):
        self.validate()
        return find_intersections(self,exclusive = True)
    
    def to_DataFrame(self):
        self.validate()
        df = pd.DataFrame(index=self.keys())
        for k,v in zip(self.keys(),self.values()):
            df.loc[k,'assigned formula'] = len(v.formula)
            df.loc[k,'mean mz'] = np.mean(v.mz)
            df.loc[k,'std mz'] = np.std(v.mz)
        return df
    
    def to_OrdinationMatrix(self, impute_value = 'nan'):
        self.validate()
        return ordination_matrix(self, impute_value = impute_value)

