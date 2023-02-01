## [1.2.3] - 01-02-2022

### Added
- msTuple.filter_bool() method

### Changed
- numpy 1.22.0 and networkx 2.5 in requirements 
- calculate_mass, aromaticity_index, element_counts, element_ratios and nominal_oxidation_state now support Cl and F (double_bond_equivalent already did)
- user_guide now lists Cl and F support and boolean indexing
- element_ratios no longer includes an optional argument for handling zero ratios (if the denominator is zero, np.nan is now returned)


## [1.2.2] - 20-09-2022

### Added
- read_csv() function
- added a to_csv() method to the msTuple class

### Changed
- msTuple.summarise() now gives min and max intensity and std mz in scientific units
- compound_class now supports elemental classification (e.g. CHO, CHON, ...)
- changed the UpSetPlot gude to be up to date with the new version of UpSetPlot.
- made small changes to the bath analysis, PCA and corems with PyKrev guides

## [1.2.1] - 04-02-2022

### Added
- average_mstuple() function

### Changed
- mass_spectrum now supports an inverted axis
- mass_spectum now plots onto a number line defined by step size
- added subset and average methods to msTupleDict 
- updated the batch analysis and user guide docs to reflect new msTupleDict and mass_spectrum behaviour
- added test for average_mstuple()
- normalise_intensity now supports row or column wise normalisation
- normalise_intensity now supports a range of transformation options, log, power2, and power3
- normalise_intensity now supports centering as a normalisation method 

## [1.2.0] - 05-01-2022 

### Added
- template function with function docstring in docs
- reaction_network and page_rank functions 
- network visualisation notebook in docs
- batch analysis notebook in docs
- msTuple and msTupleDict classes

### Changed
- all fuction docstrings edited to match the template function docstring
- all functions now take either msTuple or msTupleDict as primary input, with exception of calculate_mass which retains legacy compatibility
- aromaticity index now allowed to return negative values
- double_bond_equivalent now allowed to return negative values
- read_corems now supports P assignments 
- read_corems now checks if C,H,N,O,S,P assignments are present 
- read_corems now returns an msTuple
- git ignore now ignores windows batch commands 
- changed all relative module imports to specify the exact path
	i.e. instead of from ..formula import calculate_mass , we write from ..formula.calculate_mass import calculate_mass
- kendrick mass defect now takes string (as opposed to list) for atom_group
- no more *msTuple as a function input, replaced by msTuple_list, which is clearer and cleaner to work with
- read_batch_formularity now results in an msTuple dictionary instead of a formulaDf
- corems, upsetplot, pca and user guide notebooks to reflect 1.2.0 changes

### Deleted
- standardize_formula function as it wasn't used
- unique_formula and missing_formula, as their functionality is just a part of find_intersections and the input syntax was confusing
- unique_plot and missing_plot as they are a type of multi_van_krevelen_plot and the syntax was confusing

## [1.1.2] - 19-07-2021 

### Added
- Added __version__ to root directory init.py
- Added read_corems.py function to pykrev/formula 
- Added corems_with_pykrev user guide

### Changed
- kendrick_mass_defect doesn't require a formula list anymore 
- added two extra rounding methods for calculating nominal mass in kendrick mass defect. 
	'floor' - always round down to nearest integer, and 'ceil' always round up to nearest integer. 
- renamed the integer rounding method in kendrick mass defect from 'integer' to 'rint' to avoid confusion.
- removed requirement for formula list from kendrick mass defect plot
- made changes to tests to reflect this 
- changed readme to reflect new user guide
- changed read_formularity and read_batch_formularity to automatically exclude all formula with C13 assignments 
- changed the user guide to clarify on PyKrev's formula handling, including isotopologues and valid elements.
- changed the input and return order of filter_spectral_interference
- filter_spectral_interference now returns filter_formula as a list, not a numpy array
- monoisotopic masses in calculate_mass made more precise
- monoisotopic masses now allows you to calculate the mass of charged ions and deprotonated ions
- updated the user guide to show how mass error can be calculated using calculate mass
- read_formularity now only returns formula, intensity and mass .Note all functions should input/return in this order.
- normalise_intensity now includes pareto scaling
- addede a test case for calculate_mass to test_formularity_unittest

## [1.1.1] - 28-04-2021

### Added
- PyKrev is now published in JASMS. Added a link to the citation in read me. 

### Changed
- updated kendrick_mass_defect to support two types of rounding calculation. 
- updated tests module to test both types of kmd rounding calculation. 
- updated user guide to demonstrate both types of kmd rounding calculation. 
- updated user guide  to remove error message from bray curtis matrix section. 

## [1.1.0] - 23-03-2021

### Changed
- replaced relative_intensity with normalise_intensity, a more powerful normalisation function
- replaced code using relative_intensity in calculate diversity and mass spectrum plot 
- updated the docs to include a section on how to use PyKrev if you don't have Python installed 
- updated the user guide docs to show normalise_intensity
- updated the PCA with PyKrev guide to reflect changes to normalise_intensity
- removed all markdown documentation! Documentation is now viewed using NB viewer. This means updates to ipynb files are updated in the static docs instantly. 
- added a requirements.txt file in the root to better list dependencies (this is required by binder, but it hasn't been linked with setup.py yet.). 

## [1.0.2] -03-02-2021

### Changed
- docs organisation: each md and ipynb file is in its own folder. this makes saving the output much easier 
- pca with pykrev doc has been updated to analyse a larger dataset
- upsetplots with pykrev edited to reflect v1.0 formularity parsing 
- user guide edited with new links 

## [1.0.1] -15-01-2021

### Changed
- compound_class 'FORM' method patched so that non-assignments are appended to list and counted
- added two new loading plot examples to the pca docs 
- reformatted upsetplot docs
- short description in setup.py

## [1.0.0] -11-01-2021

### Added
- a changelog
- new plotting unittests for van krevelen patch
- read batch formularity function
- PCA with PyKrev markdown and ipynb in docs. 

### Changed
- van_krevelen_plot and multi_van_krevelen_plot now include options for compound class patches
- compound_class now includes 'FORM' alogorithm
- ordination_matrix now includes impute values option 
- default colours for multi_van_krevelen_plot now colour blind safe
- several changes to docs to showcase new functions and fixed relative intensity example
- all plotting functions now return figure and axes handles
- read formularity now gives mass error data
- read formularity now accounts for C13 assignments
- mass histogram now supports kernel density estimation overlay
- mass histogram can now plot mass error values
- typo corrections and new links in user guide

### Removed
- options for reading in peak intensity, mass, etc for read_formularity (now non optional) 
