import pandas as pd
import numpy as np 
from .msTuple import msTuple
def read_corems(corems_df, mass_type = 'calibrated', remove_multiply_assigned_peaks = True, verbose = False):
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
    mass_type: String, what form to read the mass from the dataframe in. One of: 
        'calculated' (i.e. exact formula mass),
        'calibrated' (experimental mass after calibration)
        'experimental' (raw mass from spectrum).
    remove_multiply_assigned_peaks: bool, if True remove peaks that have multiple formulae assigned to them. 
        if False msTuple will contain multiple formulae assignments
    verbose: bool, if True print number of assigned peaks and generated formulae to output

    Info
    -----------    
    Corems returns isotope assignments that pykrev cannot parse. 
    PyKrev will filter out formulae with isotopes assigned.
    """
    #Tests
    assert(mass_type in ['calculated','calibrated','experimental']), 'incorrect mass_type given'
    assert(type(remove_multiply_assigned_peaks) == bool), 'provide a boolean'
    #Setup
    assigned = corems_df['Is Isotopologue'] == 0
    assignedDf = corems_df[assigned].copy()
    total_peak_number = len(set(corems_df['Index']))
    assigned_peak_number = len(set(assignedDf['Index']))
    generated_formulae_number = assignedDf.shape[0]
    if verbose == True:
        print(f'total peaks: {total_peak_number}')
        print(f'assigned peaks: {assigned_peak_number}')
        print(f'generated formulae: {generated_formulae_number}')
        print('**************************************************')
    #deal with multiple assignments
    if generated_formulae_number > assigned_peak_number:
        if remove_multiply_assigned_peaks == True:
            didx = pd.Index([])
            for i in assignedDf['Index']:
                if len(assignedDf[assignedDf['Index'] == i]) > 1: 
                    didx = didx.union(assignedDf[assignedDf['Index'] == i].index)
            assignedDf = assignedDf.drop(index = didx)
    N = assignedDf.shape[0]
    if verbose == True and remove_multiply_assigned_peaks == True:
        print(f'{assigned_peak_number - N} multiply assigned peaks removed')
        print('----------------------------------------------------')
    #Main
    ## first check that there are columns in the dataframe corresponding to C,H,N,O,P,S,Cl and F if there aren't, we make them!
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
    if 'Cl' not in assignedDf.columns:
        assignedDf['Cl'] = [""] * N       
    if 'F' not in assignedDf.columns:
        assignedDf['F'] = [""] * N  
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
    Clstring = []
    for clnumber in assignedDf['Cl']:
        try:
            if int(clnumber) == 0:
                Clstring.append("")
            else: 
                Clstring.append(str(f"Cl{int(clnumber)}"))
        except ValueError:
            Clstring.append("")
    Fstring = [] 
    for fnumber in assignedDf['F']:
         try:
             if int(fnumber) == 0:
                 Fstring.append("")
             else: 
                 Fstring.append(str(f"F{int(fnumber)}"))
         except ValueError:
             Fstring.append("")
    pykrev_formula = [c + h + n + o + p + s + f + cl for c, h, n, o, p, s, f, cl in zip(Cstring, Hstring, Nstring, Ostring, Pstring, Sstring, Fstring, Clstring)]
    pykrev_abundance = assignedDf['Peak Height']
    if mass_type == 'calibrated':
        pykrev_mass = assignedDf['Calibrated m/z']
    elif mass_type == 'calculated':
        pykrev_mass = assignedDf['Calculated m/z']
    else:
        pykrev_mass = assignedDf['m/z']
    return msTuple(pykrev_formula, np.array(pykrev_abundance), np.array(pykrev_mass))
