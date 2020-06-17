# PyKrev: van Krevelen analysis and beyond

The aim of this package is to enable analysis of metabolomics or environmental chemistry mass-spec data in python. This includes conventional van Krevelen analysis and beyond - the package should include the most recent techniques used to analyse mass-spec data including adaptations on van-krevelen, and integrate methods used in metabolomics and environmental chemistry into one place enabling better cross-talk between the fields.

# Usage: 
See the PyKrev Examples .ipynb or .html file for examples on how the package can be used (note u must be inside the directory to import PyKrev). Additionally, each function has a help doc string that can be accessed by calling help() in python. 

Currently the package works with only two input data sets: lists of atomic formula, and lists or numpy arrays containing corresponding peak intensities. 

# Current features: 
- basic atomic formula manipulation: count atoms in formula, work out atom ratios in formula, calculate double bound equivalence, aromaticity index, and nominal oxidation state of Carbon (see: Kroll, Jesse H., et al. "Carbon oxidation state as a metric for describing the chemistry of atmospheric organic aerosol." Nature chemistry 3.2 (2011): 133.)
- comparisons of multiple lists of atomic formula to extract the unique formula in each list, and the 'missing' formula (the formula present in every other list apart from the one in question).
- calculation of various diversity indices based on formula and peak intensities. See: Mentges, Andrea, et al. "Functional molecular diversity of marine dissolved organic matter is reduced during degradation." Frontiers in Marine Science 4 (2017): 194.
- can create a variety of plots: standard van krevelen diagrams, multi_van_krevelen_diagrams (multple groups of formula one diagram), unique plots (comparing and plotting the unique formula from many lists on one diagram), missing plots (comparing and plotting the missing formula from many lists on one diagram), unique v missing subplots (subplotting the unique v missing formula from many lists on one diagram) and histogram (creating a density histogram  of the formula).

# To do:
- include compound classification function based on Rivas-Ubach, Albert, et al. "Moving beyond the van Krevelen diagram: A new stoichiometric approach for compound classification in organisms." Analytical chemistry 90.10 (2018): 6152-6160.
- create empirically dervied compound class polygons from the brite database to overlay on van Krevelen diagrams like in Brockman, Stephen A., Eric V. Roden, and Adrian D. Hegeman. "Van Krevelen diagram visualization of high resolution-mass spectrometry metabolomics data with OpenVanKrevelen." Metabolomics 14.4 (2018): 48.
- create compound class polygons based on Feunang, Yannick Djoumbou, et al. "ClassyFire: automated chemical classification with a comprehensive, computable taxonomy." Journal of cheminformatics 8.1 (2016): 61. which are more sensible than biologically derived classifications. 
- function to create upset plots and special upset plots based on compound classes, and compound class bar charts. 
- function to create upset plots that can count intersections between groups and define most significant intersections 'group similarity'
- improve the atom counts function so it can handle different styles of writing formula. 
- improve error handling in all functions by checking the validity of input variables.
- add package to python database as an importable package using https://packaging.python.org/tutorials/packaging-projects/
