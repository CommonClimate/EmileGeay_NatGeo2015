# EmileGeay_NatGeo2015
Code and data to reproduce the analysis of *Emile-Geay et al, 2015: Linkages between tropical Pacific seasonal, interannual, and orbital variability during the Holocene, Nature Geocience.
(full publications details to follow)*

Matlab code to reproduce the figures of the paper is provided in `code/matlab/``.
Running holocene_proxy_workflow at the Matlab prompt will generate .mat files from the original data files (.csv or .xls), run the block-bootstrap analysis, and produce figures. Note that the uncertainty quantification of Carré et al [2014] was carried out separately, and is stored in /data/obs/Quantiles_Carre_orig.mat and /data/obs/Quantiles_Carre_ratios.mat'.

The original data files for Holocene proxy observations are located in data/obs/.
a mat file containing the result of the block-bootstrap sampling of PMIP3 results is located in data/pmip3/.  

The relevant Python code is in code/python.

Acknowledgments go to several splendid open-source codes without which this work would have been markedly more painful, if not impossible:

**Matlab**
- [m_map](http://www.eos.ubc.ca/~rich/map.html)
- [panel](http://www.mathworks.com/matlabcentral/fileexchange/20003-panel)
- [latextable](http://www.mathworks.com/matlabcentral/fileexchange/44274-latextable)
- [wavelet toolbox](http://paos.colorado.edu/research/wavelets/)
- [Total Least Squares](http://www.mathworks.com/matlabcentral/fileexchange/31109-total-least-squares-method)
**Python**
- [seaborn](http://stanford.edu/~mwaskom/software/seaborn/)
- [anaconda](https://www.continuum.io/downloads)  

I have tried to honor licenses, but please let me know if I have strayed.
