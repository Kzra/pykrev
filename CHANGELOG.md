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