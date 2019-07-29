#' Loads a character vector of packages, installing them from Bioconductor if necessary.
#' This will also load most CRAN packages.
#'  
#' @title load_it
#' @name load_it
#' @param pack -- A vector of character strings containing package names.
#' @return nothing
#' @export
#'
load_it = function(pack, update=FALSE) {
	invisible(sapply(pack, function(package) {
		if(!require(package, quietly=T, character.only=T)){
		    if (!require("BiocManager")) install.packages("BiocManager")
		    BiocManager::install(c(package))
				library(package, quietly=T, character.only=T)
		}
	}))
} #end load_it
