# About
Fused Quadratic Inference Functions (fusQIF)

This is a repository for the R package fusQIF. The main functions divide data into blocks according to supplied indicators, fuse blocks using QIF and ADMM and re-estimate parameters using a meta-estimator asymptotically equivalent to the optimal generalized method of moments (GMM) equation. The R package's main files are:
- src/QIF-fusion.cpp: this file defines the Rcpp functions that compute the fused updates for the mean regression parameters.
- R/estimation-funcs-package.R: this file defines the R functions for the fused QIF estimation of model parameters.
- R/generate-datasets.R: this file defines the R functions to generate sample data to reproduce results in the manuscript
- the simulation_pipeline_code folder is not part of the R package but is made available to reproduce the results in Hector (2023). This folder has its own README file.

The fusQIF man file contains three examples for running the regression models from the paper.

Please email ehector@ncsu.edu with any questions or bug-reports.

# Installation

The fusQIF R package can be installed in one of two ways:
- from the downloaded gzipped tarball as R CMD INSTALL fusQIF_1.0-1.tar.gz
- from the downloaded fusQIF folder as R CMD build fusQIF and R CMD INSTALL fusQIF_1.0-1.tar.gz

Please make sure to have all packages listed in the DESCRIPTION file already installed.

# Citation

If you use the fusQIF R package, please consider citing the relevant manuscript: Emily C. Hector. Fused mean structure learning in data integration with dependence. The Canadian Journal of Statistics (2023). doi: 10.1002/cjs.11797.

# References

Emily C.Hector and Peter X.-K. Song (2020). Doubly distributed supervised learning and inference with high-dimensional correlated outcomes. Journal of Machine Learning Research, 21:1–35.
  
Emily C. Hector and Peter X.-K. Song (2020). A distributed and integrated method of moments for high-dimensional correlated data analysis. Journal of the American Statistical Association, 116:805-818.

The posdef.matrix function was written by Ravi Varadhan: https://stat.ethz.ch/pipermail/r-help/2008-February/153708.
