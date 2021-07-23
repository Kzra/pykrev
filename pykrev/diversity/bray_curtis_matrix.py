import numpy as np
from scipy import spatial
def bray_curtis_matrix(matrix):
    """ 
	Docstring for function pyKrev.bray_curtis_matrix 
	====================
	This function Scipy bray curtis function to compute a dissimilarity matrix of dimensions row * row.  
    
	Use
	----
	bray_curtis_matrix(Y)
    
	Returns a numpy array of shape(len(Y[:,0]),len(Y[:,0])) in which each value is the bray curtis dissimilarity value.  
    
	Parameters
	----------
	Y: A numpy array containing peak intensities - where rows correspond to samples and columns correspond to molecular formula
    
	Info
	----------
	The Bray-Curtis dissimilarity is always a number between 0 and 1. If 0, the two samples share all the same formula; if 1, they don’t share any formula. 
        
    """  
    assert(isinstance(matrix,np.ndarray)), 'must provide a numpy array'
    assert(len(matrix.shape) != 1), 'must provide at least two columns'
    row,col = np.shape(matrix)
    transformed_matrix = np.zeros((row,row))
    for x in range(0,row):
        for i in range(0,row):
            transformed_matrix[i,x] = spatial.distance.braycurtis(matrix[x,:],matrix[i,:])
    return transformed_matrix