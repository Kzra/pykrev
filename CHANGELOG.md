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