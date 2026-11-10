# 0. R package "etlcov"
We developed the R package 'etlcov' to run our algorithm. This package can be loaded directly in R by installing the file 'etlcov.zip'. The PDF manual for this package can be found within the file.

# 1. Overview

This repository contains the following files with their serve as below.

- **ETL-cov-simu-1**: It includes all functions used in implementing experiment 1, which is designed for general sparse covariance matrix.
- **update-simu-2**: It provides all functions used in implementing experiment 2, which is designed for bandable covariance matrix.
- **update-simu-3**: It includes all functions used in implementing experiment 3, which investigates the more challenging setting where informative sources vary across correlation entries and variance coordinates. 
- **ETL-Glasso**: This file presents an example of our method implemented for the estimation of Gaussian graphical models.
- **ETL-PCA**: This file contains the embedding functions used to implement the GB-PCA method, which serves as a baseline for comparison with the method proposed in this paper in PCA applications.
- **real data analysis**: This document provides the code and explanations for the real data analysis.

# 2. Setup Instructions

- R version: It is sufficient to ensure that the R version is R 4.3.3 or updated ones.
- Packages: All packages needed for each R program are placed at the top. All these packages can be installed using conventional methods, such as input "iinstall.packages("package name") in R.
