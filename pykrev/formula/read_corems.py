import pandas as pd
import numpy as np 
def read_corems(corems_df, mass_type = 'calibrated'):
    """ 
    Docstring for function PyKrev.read_corems
    ====================
    This function reads a pandas dataframe produced in corems using mass_spectrum_obj.to_dataframe() and outputs a list of formula, and a numpy array of peak intensities and masses.
    
    Use
    ----
    read_corems(corems_df)
    
    Returns a list of formula, a numpy array of peak intensities, a numpy array of masses.
    
    Parameters
    ----------
    mass_type: What form to read the mass from the dataframe in. 
               One of: 'calculated' (i.e. exact formula mass), 'calibrated' (experimental mass after calibration) or 'experimental' (raw mass from spectrum).
    
    Note: corems returns isotope assignments that pykrev cannot parse. 
               PyKrev will filter out formula with isotopes assigned.
    
    """
    assert(mass_type in ['calculated','calibrated','experimental']), 'incorrect mass_type given'
    assigned = corems_df['Is Isotopologue'] == 0
    assignedDf = corems_df[assigned].copy()
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
    Sstring = [] 
    for snumber in assignedDf['S']:
        try:
            if int(snumber) == 0:
                Sstring.append("")
            else: 
                Sstring.append(str(f"S{int(snumber)}"))
        except ValueError:
            Sstring.append("")
    pykrev_formula = [c + h + n + o + s for c, h, n, o, s in zip(Cstring, Hstring, Nstring, Ostring, Sstring)]
    pykrev_abundance = assignedDf['Peak Height']
    if mass_type == 'calibrated':
        pykrev_mass = assignedDf['Calibrated m/z']
    elif mass_type == 'calculated':
        pykrev_mass = assignedDf['Calculated m/z']
    else:
        pykrev_mass = assignedDf['m/z']
        
    return pykrev_formula, pykrev_abundance, pykrev_mass