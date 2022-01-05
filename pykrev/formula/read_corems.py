import pandas as pd
import numpy as np 
from .msTuple import msTuple
def read_corems(corems_df, mass_type = 'calibrated'):
    """ 
    Docstring for function PyKrev.read_corems
    ====================
    This function reads a pandas dataframe produced in corems using mass_spectrum_obj.to_dataframe() and outputs an msTuple.
    
    Use
    ----
    read_corems(corems_df)
    
    Returns an msTuple
    
    Parameters
    ----------
    mass_type: String, What form to read the mass from the dataframe in. One of: 
        'calculated' (i.e. exact formula mass),
        'calibrated' (experimental mass after calibration)
        'experimental' (raw mass from spectrum).

    Info
    -----------    
    Corems returns isotope assignments that pykrev cannot parse. 
    PyKrev will filter out formula with isotopes assigned.
    """
    #Tests
    assert(mass_type in ['calculated','calibrated','experimental']), 'incorrect mass_type given'
    #Setup
    assigned = corems_df['Is Isotopologue'] == 0
    assignedDf = corems_df[assigned].copy()
    N = assignedDf.shape[0]
    #Main
    ## first check that there are columns in the dataframe corresponding to C,H,N,O,P,S if there aren't, we make them!
    if 'C' not in assignedDf.columns:
        assignedDf['C'] = [""] * N
    if 'H' not in assignedDf.columns:
        assignedDf['H'] = [""] * N
    if 'N' not in assignedDf.columns:
        assignedDf['N'] = [""] * N      
    if 'O' not in assignedDf.columns:
        assignedDf['O'] = [""] * N           
    if 'S' not in assignedDf.columns:
        assignedDf['S'] = [""] * N       
    if 'P' not in assignedDf.columns:
        assignedDf['P'] = [""] * N      
    Cstring = []
    for cnumber in assignedDf['C']:
        try: 
            if int(cnumber) == 0:
                Cstring.append("")
            else:
                Cstring.append(str(f"C{int(cnumber)}"))
        except ValueError:
            Cstring.append("")
    Hstring = [] 
    for hnumber in assignedDf['H']:
        try: 
            if int(hnumber) == 0:
                Hstring.append("")
            else:
                Hstring.append(str(f"H{int(hnumber)}"))
        except ValueError:
            Hstring.append("")        
    Nstring = [] 
    for nnumber in assignedDf['N']:
        try:
            if int(nnumber) == 0:
                Nstring.append("")
            else:
                Nstring.append(str(f"N{int(nnumber)}"))
        except ValueError:
            Nstring.append("")
    Ostring = [] 
    for onumber in assignedDf['O']:
        try:
            if int(onumber) == 0:
                Ostring.append("")
            else:
                Ostring.append(str(f"O{int(onumber)}"))
        except ValueError:
            Ostring.append("")
    Pstring = []
    for pnumber in assignedDf['P']:
        try:
            if int(pnumber) == 0:
                Pstring.append("")
            else: 
                Pstring.append(str(f"P{int(pnumber)}"))
        except ValueError:
            Pstring.append("") 
    Sstring = [] 
    for snumber in assignedDf['S']:
        try:
            if int(snumber) == 0:
                Sstring.append("")
            else: 
                Sstring.append(str(f"S{int(snumber)}"))
        except ValueError:
            Sstring.append("")
    pykrev_formula = [c + h + n + o + p + s for c, h, n, o, p, s in zip(Cstring, Hstring, Nstring, Ostring, Pstring, Sstring)]
    pykrev_abundance = assignedDf['Peak Height']
    if mass_type == 'calibrated':
        pykrev_mass = assignedDf['Calibrated m/z']
    elif mass_type == 'calculated':
        pykrev_mass = assignedDf['Calculated m/z']
    else:
        pykrev_mass = assignedDf['m/z']
    return msTuple(pykrev_formula, np.array(pykrev_abundance), np.array(pykrev_mass))