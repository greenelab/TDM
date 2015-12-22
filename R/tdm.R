# TDM (Training Distribution Matching)
# 
# This file contains functions used to perform training distribution 
# matching on the gene expression values in a file, as well as other 
# functions that may be useful when using machine learning on such 
# datasets. It is a gentle 'normalization' technique that equalizes 
# values in the tails of the data to make the tails shorter and the 
# overall distribution more similar to those in a reference file. 
# It assumes that the target distribution is more skewed than the 
# reference distribution, such as is usually the case between RNA-seq 
# and microarray expression datasets.
# 
# author: Jeffrey Ahearn Thompson

#' Inverse Log2 of a Single Value
#' 
#' Takes the inverse log of a value.
#' 
#' @param value -- A numeric value.
#' 
#' @return the inverse log of value.
#'
#' @export
inv_log <- function(value) {
	return(2^as.numeric(value))
} # end inv_log

#' Log2 plus 1 a Single Value
#' 
#' Takes the log2 of 1 plus a value.
#' 
#' @param value -- A numeric value.
#' 
#' @return the log2(1 + value)
#'
#' @export
log2p1 <- function(value) {
	return(log1p(as.numeric(value))/log(2))
}

#' Data.table MAX
#' 
#' Finds the max of a data.table.
#' 
#' @param datatable -- A data.table whose first column will be ignored
#' but whose other columns are all numeric.
#' @return the maximum numeric value in the data.table.
#' @export
max_dt <- function(datatable) {
	return(max(as.numeric(as.matrix(datatable[,2:ncol(datatable),with=F]))))
} # end max_dt

#' Data.table MIN
#' 
#' Finds the min of a data.table.
#' 
#' @param datatable -- A data.table whose first column will be ignored
#' but whose other columns are all numeric.
#' @return the minimum numeric value in the data.table.
#' @export
min_dt <- function(datatable) {
	return(min(as.numeric(as.matrix(datatable[,2:ncol(datatable),with=F]))))
} # end min_dt

#' IQR
#' 
#' Finds the IQR of a data.table.
#' 
#' @param datatable -- A data.table whose first column will be ignored
#' but whose other columns are all numeric.
#' @return IQR of the data.table.
#' @export
iqr_dt <- function(datatable) {
	quartiles=(summary(as.numeric(data.matrix(datatable[,2:ncol(datatable), with=F])))[c(2,5)])
	return(as.numeric(quartiles[2]) - as.numeric(quartiles[1]))
} # end iqr_dt

#' Quartiles
#' 
#' Finds the 1st and 3rd quartiles of a data.table.
#' 
#' @param datatable -- A data.table whose first column will be ignored
#' but whose other columns are all numeric.
#' @return 1st and 3rd quartiles of the data.table.
#' @export
quartiles_dt <- function(datatable) {
	quartiles=(summary(as.numeric(array(as.matrix(datatable[,2:ncol(datatable), with=F]))))[c(2,5)])
	return(quartiles)
}

#' Inverse Log Transform
#' 
#' Given a data.table of values, each value is inverse log transformed.
#' 
#' @param data -- A data.table of values. The first column should be
#' gene symbols and the rest should be expression values.
#' @param file -- A filename with expression values. This is an alternative
#' to specifying data.
#' @return the transformed data.table.
#' 
#' @export
inv_log_transform = function(data = NULL, file = NULL) {
	if(is.null(data)) {
		if(is.null(file)) {
			stop("You must specifcy either data or file.")
		} else {
			data <- fread(file, header=T, data.table=T)
		}
	} 
	
	# Convert the numeric part of the data.table to a matrix.
	# Assuming first column contains gene symbols.
	datamat = data.matrix(data[,2:ncol(data),with=F])
	
	# Inverse log transform the matrix.
	inv_log = apply(datamat,1,inv_log)
	
	# Convert the result back to a data.table and bind the symbols back on.
	result = data.table(cbind(data[[1]], t(inv_log)))
	setnames(result, colnames(result), colnames(data))
	
	return(result)
} # end inv_log_transform

#' Log2 Transform 1 Plus
#' 
#' The input is log transformed such that all values are non-negative by
#' adding 1 to each value before transformation.
#' 
#' @param data -- A data.table of values. The first column should be
#' gene symbols and the rest should be expression values.
#' @param file -- A filename with expression values. This is an alternative
#' to specifying data.
#' 
#' @return the transformed data.table.
#' 
#' @export
log_transform_p1 = function(data = NULL, file = NULL) {
	if(is.null(data)) {
		if(is.null(file)) {
			stop("You must specify either data or file.")
		} else {
			data <- fread(file, header=T, data.table=T)
		}
	} 
	
	# Convert the numeric part of the data.table to a matrix.
	# Assuming first column contains gene symbols.
	datamat = data.matrix(data[,2:ncol(data),with=F])
	
	# Log (+1) transform the matrix.
	log_transform = apply(datamat,1,log2p1)
	
	# Convert the result back to a data.table and bind the
	# gene symbols back on.
	result = data.table(cbind(as.character(data[[1]]), t(log_transform)))
	setnames(result, colnames(result), colnames(data))
	
	return(result)
} # end log_transform_p1

#' Zero to One Transform a Vector
#' 
#' A data vector is transformed so that all values are in range [0,1].
#' 
#' @param data -- A vector of gene expression values.
#' 
#' @return the transformed vector.
#' 
#' @export
zero_to_one = function(data) {
	# If the data are already all 0, then do nothing.
	# If the data are all the same, then stop.
	if(sum(data)==0) {
		return(as.vector(rep(0.0, length(data))))
	}else if(min(data) == max(data)){
		return(as.vector(rep(0.0, length(data))))
		if(min(data) < 1){
			return(data)
		} else {
			return(as.vector(rep(0.0, length(data))))
		}
	}
	
	# Transform the values.
	line = (data-min(data))/(max(data)-min(data))
	
	return(as.vector(line))
}

#' Zero to One Transformation
#' 
#' The input is standardized to the range [0,1] by row.
#' 
#' @param data -- A data.table of values. The first column should be
#' gene symbols and the rest should be expression values.
#'  
#' @return the transformed datatable.
#' 
#' @export
zero_to_one_transform = function(datatable) {
	# Convert the data to a matrix.
	datamat = data.matrix(datatable[,2:ncol(datatable),with=F])
	
	# Transform to [0,1] range, row by row.
	zo = apply(datamat,1,zero_to_one)
	
	# Bind on the gene symbols.
	result = data.table(data.frame(datatable[[1]], t(zo)))
	setnames(result, colnames(result), colnames(datatable))
	
	return(result)
}

#' TDM Transformation on a Single Value
#' 
#' Applies TDM transformation to a single value.
#' 
#' @param value -- A single value from a distribution that is being
#' TDM transformed.
#' @param iqr -- inter-quartile-range of distribution being transformed
#' @param old_min -- min of distribution being transformed
#' @param old_max -- max of distribution being transformed
#' @param old_first_q -- first quartile of distribution being transformed
#' @param old_third_q -- third quartile of distribution being transformed
#' @param ref_min -- min of reference distribution
#' @param ref_max -- max of reference distribution
#' @param ref_third_q -- third quartile of reference distribution
#' @param ref_first_q -- first quartile of reference distribution
#' @param log_reference -- logical indcating if reference distribution was logged
#' before values were extracted
#' @param log_target -- logical indicating if target distribution should be logged
#' 
#' @return TDM of the value passed.
#' @export
tdm <- function(value, iqr, old_min, old_max, old_third_q, old_first_q, ref_min, ref_max, ref_third_q, ref_first_q, inv_reference=TRUE, log_target=TRUE) {

	# If the reference distribution was log transformed, then inverse log its
	# values here. This will come out the same, because log is monotonically
	# increasing.
	if(inv_reference) {
		ref_max = 2^ref_max
		ref_min = 2^ref_min
		ref_third_q = 2^ref_third_q
		ref_first_q = 2^ref_first_q
	}
	
	# Find the number of iqrs from reference distribution that fit between
	# it's third quartile and max.
	topscale = (ref_max- ref_third_q)/(ref_third_q - ref_first_q)
	
	# Do the same for min and first quartile.
	bottomscale = (ref_first_q - ref_min)/(ref_third_q - ref_first_q)
	
	# Set cutoff values that would make the target distribution have same
	# relationship between quartiles and min and max.
	upandout = old_third_q + topscale*iqr
	downandout = old_first_q - bottomscale*iqr
	
	# The bottom threshold should not be less than the minimum already was.
	if(downandout < old_min) {
		downandout = old_min
	}

	# Set any values that are greater or less than the thresholds to be equal
	# to the thresholds.
	if(value > upandout) {
		value = upandout
	} else if(value < downandout) {
		value = downandout
	}
	
	# Now scale the value so that all values are within the same min and max
	# as the reference distribution.
	result = (value - downandout)/(upandout - downandout) * (ref_max - ref_min) + ref_min
	
	# If the target distribution should be log transformed, then do it.
	if(log_target) {
		result = log2(result)
	}
	
	return(result)

} # end tdm

euc.dist = function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

trim = function (x) gsub("^\\s+|\\s+$", "", x)

#' TDM Transformation
#' 
#' The expression files should have the same sorts of gene symbols. The target file
#' will be filtered to include only symbols that are in the reference file and the
#' order will also be changed to match. If a gene is not present in the target file,
#' an entry of all 0s will be added for that gene.
#' 
#' @param file -- The file of gene expression values to transform, presumably
#' containing RNA-seq values. It should have a header, but the header
#' should only contain column names, not a row name. All other rows should
#' have the same number of values and the first value should be a row name.
#' This is an inconveniance of working with data.table, but they are much faster.
#' Values should be tab-separated.
#' 
#' @param ref_file -- the file to transform in relation to, presumably microarray
#' expression values. Values should be tab-separated.
#' 
#' @param negative -- should the reference file be first inverse log transformed, then
#' log transformed again adding 1 to each value? Do this if the file has already
#' been log transformed and their are negative numbers in it.
#' 
#' @param filter_p -- should the patients in the target file be filtered to include
#' only those that match patients (column ids) in the reference file? Only the first
#' 15 characters are checked for match.
#' 
#' @param inv_reference -- should reference data be inverse log transformed? Do this 
#' if the target is not already log transformed (which is the typical case). This just
#' ensures that the comparison is meaningful.
#' 
#' @param log_target -- should the target file be log transformed after TDM? Usually
#' it should, since microarray data are typically log transformed and the whole
#' point of the transformation is to make the data more comparable.
#' 
#' @return void
#'
#' @export
tdm_transform <- function(target_data=NULL, 
	ref_data=NULL, 
	file=NULL, 
	ref_file=NULL, 
	negative=FALSE, 
	filter_p=FALSE, 
	inv_reference = TRUE, 
	log_target=TRUE){
	
	load_it(c("data.table", "binr", "scales"))

	# Preprocess the reference data.
	if(is.null(ref_data)) {
		# Read the first line of the reference file.
		ref_head = readLines(ref_file, n=1)
		
		# Split all values in the header of the reference file by tabs.
		split_head = unlist(strsplit(ref_head, "\t"))
	
		# Read in the reference expression values.
		ref_values <- fread(ref_file, header=F, skip=1, data.table=T)
		
		# If there is no row name for the header, then handle it.
		if(ncol(ref_values) == length(split_head)) {
			ref_head = split_head[2:length(split_head)]
		} else {
			ref_head = split_head
		}
		rm(split_head)
		
		# Add a rowname to the header, just to make things consistent.
		setnames(ref_values, colnames(ref_values), c("gene", ref_head))
		rm(ref_head)
	} else {
		ref_values = ref_data
	}

    # Get the gene symbols from the reference data.
	genes = data.frame(gene = ref_values[[1]], drop=F)
	
	# Inverse log, then relog reference values if asked.
	if(negative) {
		ref_values = inv_log_transform(ref_values)
		ref_values = log_transform_p1(ref_values)
		ref_values$gene = genes$gene
	}
	
	# Key the table to gene symbol.
	#setkey(ref_values, "gene")

	# Preprocess the target data.
	if(is.null(target_data)) {
		# Get the header of the target file.
		exp_head = readLines(file, n=1)
		split_head = unlist(strsplit(exp_head, "\t"))
		
		# Read in the target expression file.
		expression_values <- fread(file, header=F, skip=1, data.table=T)
		
		# Again, handle the case when there is no row name for the target
		# file header.
		if(ncol(expression_values) == length(split_head)) {
			exp_head = split_head[2:length(split_head)]
		} else {
			exp_head = split_head
		}
		rm(split_head)
		setnames(expression_values, colnames(expression_values), c("gene", exp_head))
		rm(exp_head)
	}else {
		expression_values = target_data
	}
	
	# Key the table to gene symbol.
	#setkey(expression_values, "gene")
	
	# If the gene symbols come with a vertical bar and Entrez id,
	# such as Illumina data often do, then strip that suffix.
	expression_values$gene = gsub('^(.*)\\|.*$', '\\1', expression_values$gene)
	
	expression_values$gene = trim(expression_values$gene)

	# Figure out which genes are in the reference file, but not
	# in the target file.
	missing = ref_values$gene[!(ref_values$gene %in% expression_values$gene)]
	
	# Order the genes in the target file to be the same as those in
	# the reference file.
	expression_values = expression_values[match(ref_values$gene, expression_values$gene),1:ncol(expression_values),with=FALSE]

	if(nrow(expression_values) < 1) {
		stop("No matching genes found between datasets.")
	}
	
	# Filter the genes in the target file to include only those in
	# the reference file.
	expression_values = expression_values[expression_values$gene %in% ref_values$gene,1:ncol(expression_values),with=FALSE]
	
	# Add missing data entries of all 0's for the genes that were in the 
	# reference file but not in the target file.
	missing_matrix = matrix(rep(0, 
			(ncol(expression_values)-1) * length(missing)), 
			nrow=length(missing))
	
	if(nrow(missing_matrix) > 0) {
		missing_dt <- data.table(cbind(gene = missing, data.frame(missing_matrix)))
		setnames(missing_dt, colnames(missing_dt), colnames(expression_values))
		expression_values <- rbindlist(list(expression_values, missing_dt))
	}
	
	if(filter_p) {
		# Filter the patients in the target file to include only those in
		# the reference file. Currently, this assumes TCGA style patient ids.
		setnames(expression_values, colnames(expression_values), substr(colnames(expression_values),1,15))
		setnames(ref_values, colnames(ref_values), substr(colnames(ref_values),1,15))
		cols = intersect(colnames(expression_values), colnames(ref_values))
		expression_values <- expression_values[,na.omit(cols),with=FALSE]
	}

	# Set variables that are needed for TDM.
	old_min = as.numeric(min_dt(expression_values))
	old_max = as.numeric(max_dt(expression_values))
	
	old_quartiles <- quartiles_dt(expression_values)
	old_third_q <- as.numeric(old_quartiles[2])
	old_first_q <- as.numeric(old_quartiles[1])
	
	iqr <- as.numeric(iqr_dt(expression_values))
	
	ref_min <- as.numeric(min_dt(ref_values))
	ref_max <- as.numeric(max_dt(ref_values))
	
	ref_quartiles <- quartiles_dt(ref_values)
	ref_third_q <- as.numeric(ref_quartiles[2])
	ref_first_q <- as.numeric(ref_quartiles[1])

	# Get the coloumn names of the target file.
	cols <- colnames(expression_values[,2:ncol(expression_values), with=F])

	# Convert all expression entries to numeric.
	for(j in cols) set(expression_values, j=j, value=as.numeric(expression_values[[j]]))
	
	# Perform TDM
	datamat = data.matrix(expression_values[,2:ncol(expression_values),with=F])
	tdmmat = apply(datamat,1,function(x) {sapply(x, function(y) tdm(y, iqr, old_min, old_max, old_third_q, old_first_q, ref_min, ref_max, ref_third_q, ref_first_q, inv_reference, log_target))})
	result = data.table(cbind(as.character(expression_values[[1]]), t(tdmmat)))
	result[[1]] = as.factor(result[[1]])
	setnames(result, colnames(result), colnames(expression_values))
	expression_values = result

	# Order the genes in the target file to be the same as those in
	# the reference file.
	expression_values <- expression_values[match(genes$gene, expression_values$gene),1:ncol(expression_values),with=FALSE]

	return(expression_values)
} # end tdm_transform

