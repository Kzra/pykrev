import numpy as np
import pandas as pd
def relative_intensity(matrix):
    
    
    """ 
	Docstring for function pyKrev.relative_intensity 
	====================
	This function takes a matrix of peak intensities and calculates the relative intensities for each column. 
    
	Use
	----
	relative_intensity(Y)
    
	Returns a numpy array of shape(Y) in which each value is the relative intensity of that column.  
    
	Parameters
	----------
	Y: A list or numpy array containing peak intensities.  
    
    
	Info
	----------
	Relative peak intensity is calculated by summing all the intensities given in a column of Y,
	and dividing each item in the column (i) by the sum. So that sum(i...n) = 1.  
        
    """  
    
    if isinstance(matrix,pd.DataFrame): #if the user has supplied a dataframe e.g. produced by ordination_matrix
        matrix = matrix.to_numpy()
    
    reshape = False
    if len(matrix.shape) == 1: #i.e. if the matrix is 1D
        matrix = np.reshape(matrix,(len(matrix),1)) #reshape to 2D 
        reshape = True 
    row, col = np.shape(matrix)
    transformed_matrix = np.zeros((row,col))
    for i in range(0,col):
        col_sum = sum(matrix[:,i])
        for x in range(0,row):
                transformed_matrix[x,i] = matrix[x,i] / col_sum
    if reshape == True: #if matrix provided as 1D
        transformed_matrix = np.reshape(transformed_matrix,(len(transformed_matrix),)) #reshape back to 1D
    return transformed_matrix