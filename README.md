# TDM

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.32852.svg)](http://dx.doi.org/10.5281/zenodo.32852)

An R package for normalizing RNA-seq data to make them comparable to microarray data for use with machine learning.

In order to include this package in your R code, you only need to do the following:

install.packages("devtools")  
library(devtools)  
devtools::install_github("greenelab/TDM")  
library(TDM)

Alternatively, if you have latex installed, you can build the vignette:

devtools::install_github("greenelab/TDM", build_vignettes = TRUE)

A copy of the compiled vignette is also include in this repository.

Acknowledgements: 
This research is funded in part by the Gordon and Betty Moore Foundation’s Data-Driven Discovery Initiative through Grant GBMF4552 to CSG. JT is a Neukom Graduate Fellow supported by the William H. Neukom 1964 Institute for Computational Science. This work was supported in part by P20 GM103534, P30 CA023108 and UL1 TR001086 from the NIH and an American Cancer Society Research Grant, #IRG-82-003-27. The funders had no role in study design, data collection and analysis, decision to publish, or preparation of the manuscript.
