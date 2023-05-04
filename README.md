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

## FAQ

**How should the RNA-seq input data be prepared?**

*TDM expects RNA-seq count data by default. In the tutorial, the reason the counts are not discrete values is that they are expected counts from RSEM. This is also fine to use. We have not tested TDM with RPKM or other normalized values.*

**Can I use TPM corrected values with TDM?**

*Unfortunately, this is unlikely to work well. We assumed raw or RSEM counts with the TDM approach.*

**Can I use TDM to normalize array data to match RNA-seq data?**

*We generally do not advise this study design. We expect array data to have less precision at higher expression levels due to saturation, while counts-based RNA-seq data does not have that problem. We recommend reshaping the data expected to have more dynamic range (RNA-seq) to fit the narrower and less precise (array) distribution.*

