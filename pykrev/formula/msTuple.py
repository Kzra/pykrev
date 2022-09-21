from typing import NamedTuple
import numpy as np
import pandas as pd
class msTuple(NamedTuple):
    """ 
    Docstring for class pykrev.msTuple
    ==========
    Core DataType for PyKrev. A tuple representing an assigned mass spectrum and contains three objects:
            1. msTuple.formula: a list of molecular formula strings 
            2. msTuple.intensity: a numpy.ndarray of peak intensities
            3. msTuple.mz a numpy.nd array of calibrated mz values
    
    Use
    ----------
    msTuple(Y,Z,X) 
    
    Returns an msTuple
    
    Parameters 
    ----------
    Y: A list of molecular formula strings. 
    Z: A numpy.ndarray of peak intensities
    X: A numpy.ndarray of mz values.  
    
    Methods 
    ----------
    msTuple.validate(): validate the type and length of objects in an msTuple
    
    msTuple.summary(): summarise the objects in the tuple

    msTuple.filter_mz(low,high): returns a new msTuple filtered between low and high mz

    msTuple.filter_intensity(low,high): returns a new msTuple filtered between low and high intensity

    msTuple.filter_spectral_interference(): returns a new msTuple which has been filtered of spectral interference by doubley charged molecular ions (see pykrev.filter_spectral_interference)

    msTuple.to_csv(): writes the msTuple to a .csv file
    """
    
    formula: list
    intensity: np.ndarray
    mz: np.ndarray
    
    def __repr__(self) -> str:
        return f'msTuple(formula={np.array(self.formula)}), intensity={self.intensity}, mz={self.mz}'
    
    def validate(self):
        assert type(self.formula) == list, "msTuple.formula must be provided as a list"
        assert type(self.intensity) == np.ndarray, "msTuple.intensity must be provided as an ndarray"
        assert type(self.mz) == np.ndarray, "msTuple.mz must be provided as an ndarray"
        assert len(self.formula) > 1, "msTuple.formula must not be empty"
        assert len(self.formula) == len(self.intensity) == len(self.mz), "msTuple.formula, msTuple.intensity and msTuple.mz must be of equal length"
    
    def summary(self):
        self.validate()
        print(f'assigned formula = {len(self.formula)} ')
        print(f'min intensity = {np.min(self.intensity):.1E}')
        print(f'max intensity = {np.max(self.intensity):.1E}')
        print(f'mean mz = {np.mean(self.mz)} ')
        print(f'std mz = {np.std(self.mz):.1E}')
    
    def filter_spectral_interference(self, tol = 2, verbose = True):
        #Setup
        self.validate()
        mass_list = self.mz
        formula_list = np.array(self.formula)
        peak_intensities = self.intensity
        mass_defect = mass_list - np.floor(mass_list)
        c13d2 = (13.00335-12)/2
        flags = np.zeros(len(mass_list))
        #Main
        for i in range(0,len(mass_list)):
            if mass_defect[i] > 0.4 and mass_defect[i] < 0.8:
                flags[i] = 1
                mono_mass = mass_list[i] - c13d2
                mass_diff = abs(mass_list - mono_mass)
                mono_pos = mass_diff == min(mass_diff)
                if min(mass_diff)/mass_list[i] * 1e6 < tol:
                    if peak_intensities[mono_pos]/peak_intensities[i] < 10:
                        flags[mono_pos] = 2
        spectralFilter = flags == 0 
        filtermz = mass_list[spectralFilter]
        filterformula = formula_list[spectralFilter]
        filterintensity = peak_intensities[spectralFilter]
        if verbose == True:
            print(f"{len(mass_list)-len(filtermz)} interferences removed.")
        filterformula = list(filterformula) #reconvert filter_formula to list to sort out data type issue
        return self._replace(formula = filterformula, intensity = filterintensity, mz = filtermz)
    
    def filter_mz(self,low,high):
        self.validate()
        highpassBool = self.mz > low
        lowpassBool = self.mz < high
        bandpassBool = highpassBool & lowpassBool
        fArray = np.array(self.formula)
        filtermz = self.mz[bandpassBool]
        filterintensity = self.intensity[bandpassBool]
        filterformula = list(fArray[bandpassBool])
        return self._replace(formula = filterformula, intensity = filterintensity, mz = filtermz)
    
    def filter_intensity(self,low,high):
        self.validate()
        highpassBool = self.intensity > low
        lowpassBool = self.intensity < high
        bandpassBool = highpassBool & lowpassBool
        fArray = np.array(self.formula)
        filtermz = self.mz[bandpassBool]
        filterintensity = self.intensity[bandpassBool]
        filterformula = list(fArray[bandpassBool])
        return self._replace(formula = filterformula, intensity = filterintensity, mz = filtermz)

    def to_csv(self, path):
        self.validate()
        csv = pd.DataFrame()
        csv['formula'] = self.formula
        csv['intensity'] = self.intensity
        csv['m/z'] = self.mz
        csv.to_csv(path, index= False)
        print(f'msTuple written to {path}')