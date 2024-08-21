# fusQIF
Fused Quadratic Inference Functions

This is a repository for the R package to perform regression analysis of high-dimensional correlated responses. The approach divides data into blocks according to the supplied indicator, analyses blocks using quadratic inference functions, and then fuses blocks using ADMM and re-estimates parameters using optimal generalized method of moments (GMM) equation. The R package's main files are:
- src/QIF-fusion.cpp: this file defines the Rcpp functions that compute the objective functions.
- R/generate-datasets.R: this file defines the user-facing R functions for generating simulated data.
- R/estimation-funcs-package.R: this file defines the user-facing R function for the divide-and-conquer estimation of model parameters.

The fusQIF man file contains an example for running the regression model from the paper.

Please email ehector@ncsu.edu with any questions or bug-reports.

# Installation

The fusQIF R package can be installed in one of two ways:
- from the downloaded gzipped tarball as R CMD INSTALL fusQIF_1.0-1.tar.gz
- from the downloaded and renamed fusQIF folder as R CMD build fusQIF and R CMD INSTALL fusQIF_1.0-1.tar.gz

Please make sure to have all packages listed in the DESCRIPTION file already installed. If you encounter a library not found error for lgfortran, please try installing gfortran from here: https://cran.r-project.org/bin/macosx/tools/.

# Citation

If you use the fusQIF R package, please consider citing the relevant manuscript: Emily C. Hector. Fused mean structure learning in data integration with dependence. The Canadian Journal of Statistics (2023). doi: 10.1002/cjs.11797.

# References

Emily C. Hector and Peter X.-K. Song (2020). Doubly distributed supervised learning and inference with high-dimensional correlated outcomes. Journal of Machine Learning Research, 21:1–35.
  
Emily C. Hector and Peter X.-K. Song (2020). A distributed and integrated method of moments for high-dimensional correlated data analysis. Journal of the American Statistical Association, 116:805-818.
