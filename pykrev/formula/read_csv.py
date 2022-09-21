import pandas as pd
import numpy as np
from .msTuple import msTuple
def read_csv(Y, column_headers = False):
    """ 
    Docstring for function pykrev.read_csv
    ==========
    Reads an assigned mass list from a .csv file into an msTuple. 
    
    Use 
    ----------
    read_csv(Y)
    
    Returns an msTuple. 
    
    Parameters 
    ----------
    Y:  String, path to the .csv file. 
        The .csv file should have three columns in the following order: formula, intensity and mass.
    
    column_headers: Bool, does the .csv file include column headers?
    """
    if column_headers == False:
        data = pd.read_csv(Y, header = None)
    else:
        data = pd.read_csv(Y)
    msTupleObj = msTuple(list(data.iloc[:,0]),np.array(data.iloc[:,1]),np.array(data.iloc[:,2]))
    return msTupleObj