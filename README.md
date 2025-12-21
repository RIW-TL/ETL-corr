# 1. Overview

This repository contains the following files with their serve as below.

- **Main-functions.R**: It includes all functions used in our paper, such as the data generation, ETL-corr, ETL-cov and the transfer-enhanced variance estimation. The specific role of each function can be found in the corresponding inline anotations.
- **Example of ETL-corr & ETL-cov.R**: It provides an example of running the above functions and specifically how to use our method. In addition, this example contains the other competitors.
- **Functions of Trans-CLIME.R**: This file contains the embedding functions used to implement the Trans-CLIME method, which serves as a baseline for comparison with the method proposed in this paper in GGM applications.
- **Example of ETL-Glasso.R**: This file presents an example of our method implemented for the estimation of Gaussian graphical models.
- **Functions of GB-PCA**: This file contains the embedding functions used to implement the GB-PCA method, which serves as a baseline for comparison with the method proposed in this paper in PCA applications.
- **Example of ETL-PCA.R**: This file presents an example of our method implemented in the application of PCA.
- **real data analysis**: This document provides the code and explanations for the real data analysis.

# 2. Setup Instructions

- R version: It is sufficient to ensure that the R version is R 4.3.3 or updated ones.
- Packages: All packages needed for each R program are placed at the top. All these packages can be installed using conventional methods, such as input "iinstall.packages("package name") in R.

# 3. Data Acquisition

The stock data set, referred to as "stockdata" in our real-world data analysis, is available in the R package "huge". Users can access the dataset via the following methods:

library(huge)  
data("stockdata")
