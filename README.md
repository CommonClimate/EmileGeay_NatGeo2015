# EmileGeay_NatGeo2015
Code and data to reproduce the analysis of Emile-Geay et al, 2015: Linkages between tropical Pacific seasonal, interannual, and orbital variability during the Holocene, Nature Geocience.
(full publications details to follow)

Matlab code to reproduce the figures of the paper is provided in code/matlab/.
Running holocene_proxy_workflow at the Matlab prompt will generate .mat files from the original data files (.csv or .xls), run the block-bootstrap analysis, and produce figures. Note that the uncertainty quantification of Carré et al [2014] was carried out separately, and is stored in /data/obs/Quantiles_Carre_orig.mat and /data/obs/Quantiles_Carre_ratios.mat'.

The original data files for Holocene proxy observations are located in data/obs/.
a mat file containing the result of the block-bootstrap sampling of PMIP3 results is located in data/pmip3/.  

The Python code that generated it is in code/python.
