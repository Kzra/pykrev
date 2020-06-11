def relative_intensity(intensity_list):
    
    
    """ 
	Docstring for function pyKrev.relative_intensity 
	====================
	This function takes a list of peak intensities and calculates the relative intensities. 
    
	Use
	----
	relative_intensity(Y)
    
	Returns a list of len(Y) in which each item is the relative intensity.  
    
	Parameters
	----------
	Y: A list or numpy array containing peak intensities.  
    
    
	Info
	----------
	Relative peak intensity is calculated by summing all the intensities given in Y,
	and dividing each item (i) by the sum. So that sum(i...n) = 1.  
	 
        
    """  
    
    rel_intensity= []
    sum_intensity = sum(intensity_list)
    for i in intensity_list:
        rel_intensity.append(i/sum_intensity)
    return rel_intensity
