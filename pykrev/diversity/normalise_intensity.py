import numpy as np
import pandas as pd
def normalise_intensity(input_matrix, norm_method = 'sum', norm_subset =  'ALL', p_L = 500, p_P = 0.5, log = False):
    
    
    """ 
	Docstring for function pyKrev.normalise_intensity 
	====================
	This function takes an intensity data matrix and applies a normalisation method on the rows of the data.
	Normalisation consists of two processes: (1) applying a subset method to the data, (2) generating normalisation factors (e.g. mean, median, zscore)
	based on that subset which are then applied to the entire dataset.

	Use
	----
	normalise_intensity(Y)
    
	Returns a numpy array or pd.dataframe of shape(Y) in which each value corresponds to the row normalised intensity.  
    
	Parameters
	----------
	Y: A numpy array of shape Y[samples,formula] containing peak intensities OR an ordination matrix produced by pk.ordination matrix 
	method: A string decribing the relative intensity metric to be used. One of:
            - 'zscore': zscore normalisation i.e. mean(Y[i,:]) == 0, sd(Y[i,:]) == 1
            - 'pareto': pareto normalisation i.e. mean(Y[i,:]) == 0 [similar to z-score, but divide by the square root of the sd]
            - 'minmax': min max normalisation i.e. min(Y[i,:] == 0, max(Y[i,:] == 1
            - 'mean': mean normalization i.e. mean(Y[i,:]) == 0, min(Y[i,:]) >= -1, max(Y[i,:]) <= 1
            - 'median': median normalisation i.e. median(Y[i,:]) == 0, min(Y[i,:]) >= -1, max(Y[i,:]) <= 1
            - 'sum': sum relative intensity (sum(Y[i,:]) == 1)
            - 'max': max relative intensity (max(Y[i,:]) == 1)
            - 'unit_vector': unit vector sum relative intensity (sum(Y[i,:].^2) == 1)
            - 'binary': create a presence/absence matrix of formula, ignoring peak intensity.
            - 'none': do not do normalise, return the original, or log(original) input matrix.
	norm_subset: A string decribing the subset method to be used to generate the normalisation factors (tgnf). One of:
            - 'ALL': Use all formula in the dataset tgnf. 
            - 'LOS': take the top L order formula (ranked by intensity) from each sample tgnf, where p_L is the number of formula to take (default 500). 
            - 'PPP': take the proportion of peaks with a minimum percentage of observed values tgnf, where p_P defines the minimum percentage (default 0.5 (50%)). 
    p_L: parameter L in LOS subset, i.e. the number of top formula to take from each sample
    p_P: parameter P in PPP subset, i.e. the minimum percentage of observed intensities required to maintain a formula
    log: boolean, if true perform log transformation on peak intensities before subset and normalisation. if true, input_matrix should not contain zero values.
        
	Info
	----------
	For a description of unit vector sum relative intensity see:
	"On the Normalization of a Mass Spectrum for Comparison of Two Spectra - Zeev B. Alfassi ASMS (2004)"
        For a description of the other normalisation approaches used in this function, please refer to :
	"FT-ICR-MS Peak Intensity Normalization for Complex Mixture Analyses" Thompson et al. 2021.
        and "pmartR: Quality Control and Statistics for Mass Spectrometry-Based Biological Data", Stratton et al. 2019"
    """
    
    ## CHECK THE INPUT METHODS
    assert norm_method in ['zscore','minmax','mean','median','sum','max','unit_vector','binary','none', 'pareto'], "method not recognised"
    assert norm_subset in ['ALL','LOS','PPP'], "subset method not recognised"

    ## TRANSFORM THE INPUT DATA
    ordination_supplied = False
    if isinstance(input_matrix,pd.DataFrame): #if the user has supplied a dataframe e.g. produced by ordination_matrix
        ordination_supplied = True #later we will transform the data back to a dataframe
        ordination_copy = input_matrix.copy()
        input_matrix = input_matrix.to_numpy(dtype = float)
    oneDreshape = False
    if len(input_matrix.shape) == 1: #i.e. if the matrix is 1D
        input_matrix = np.reshape(input_matrix,(1,len(input_matrix))) #reshape to 2D 
        oneDreshape = True #later we will transform the data back to 1D
    if log == True: #log transform the data
        assert 0 not in input_matrix, "log 0 undefined, consider imputing 0 values as 1"
        input_matrix = np.log(input_matrix)
    
    ##CREATE THE TRANSFORMED MATRIX 
    rows,cols = np.shape(input_matrix)
    transformed_matrix = np.zeros((rows,cols))

    ## POPULATE THE TRANSFORMED MATRIX
    if norm_method == 'none':
        transformed_matrix = input_matrix
    else:
        if norm_subset == 'PPP':
            # first calculate the columns above the PPP threshold
            threshold = int(p_P * rows)
            # create a booean array corresponding to the above p_P cols in each row
            boolean = np.array([False] * cols)
            for c in range(0,cols):
                if sum(input_matrix[:,c] > 0) > threshold:
                    boolean[c] = True
            assert sum(boolean) > 0, 'you have set the p_P value to high, there are no data points in your subset'
            print(f" There are {sum(boolean)} peaks/formula in your subset.")
        else:
            # create a boolean array of all columns 
            boolean = np.array([True] * cols)
        for i in range(0,rows):
            # create a boolean array corresponding to the top p_L values in this row
            if norm_subset == 'LOS':
                dsc_idx = np.argsort(input_matrix[i,:])[::-1] #sort the values in descending order
                dsc_idx = np.roll(dsc_idx,-np.count_nonzero(np.isnan(input_matrix[i,:]))) #put nans to the end of the array
                L_order = dsc_idx[0:p_L]
                boolean = np.array([False] * cols)
                for L in L_order:
                    boolean[L] = True #populate boolean array with top L values
            # calculate the appropriate normalisation factors                         
            row_sum = np.nansum(input_matrix[i,boolean])
            row_max = np.nanmax(input_matrix[i,boolean])
            row_min = np.nanmin(input_matrix[i,boolean])
            row_mean = np.nanmean(input_matrix[i,boolean])
            row_median = np.nanmedian(input_matrix[i,boolean])
            row_std = np.nanstd(input_matrix[i,boolean])
            row_norm = np.sqrt(np.nansum(input_matrix[i,boolean]**2))
            if norm_method == 'sum':
                for x in range(0,cols):
                        transformed_matrix[i,x] = input_matrix[i,x] / row_sum
            elif norm_method == 'max':
                for x in range(0,cols):
                        transformed_matrix[i,x] = input_matrix[i,x] / row_max
            elif norm_method == 'unit_vector':
                  for x in range(0,cols):
                        transformed_matrix[i,x] = input_matrix[i,x] / row_norm
            elif norm_method == 'zscore':
                  for x in range(0,cols):
                        transformed_matrix[i,x] = (input_matrix[i,x] - row_mean)/ row_std
            elif norm_method == 'pareto':
                  for x in range(0,cols):
                        transformed_matrix[i,x] = (input_matrix[i,x] - row_mean)/ np.sqrt(row_std)
            elif norm_method == 'minmax':
                  for x in range(0,cols):
                        transformed_matrix[i,x] = (input_matrix[i,x] - row_min)/ (row_max - row_min)
            elif norm_method == 'mean':
                  for x in range(0,cols):
                        transformed_matrix[i,x] = (input_matrix[i,x] - row_mean)/ (row_max - row_min)
            elif norm_method == 'median':
                  for x in range(0,cols):
                        transformed_matrix[i,x] = (input_matrix[i,x] - row_median)/ (row_max - row_min)
            elif norm_method == 'binary':
                    row_binary = input_matrix[i,:] > 0 
                    transformed_matrix[i,row_binary] = 1

    ## RETRANSFORM THE DATA AND RETURN
    if oneDreshape == True: #if input_matrix provided as 1D
        transformed_matrix = np.reshape(transformed_matrix,(transformed_matrix.shape[1],)) #reshape back to 1D
    if ordination_supplied == True: #if input_matrix supplied as pk.ordination_matrix
        ordination_copy.iloc[:,:] = transformed_matrix #repopulate the ordination matrix with transformed values
        return ordination_copy #return an ordination matrix
    else:
        return transformed_matrix #else return a numpy array
