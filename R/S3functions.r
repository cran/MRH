###########################################################################################
# Contains the functions needed to: convert to an MRH object, summary, and plot for PH model.
###########################################################################################
library(survival)
library(coda)

MRH = function(x) UseMethod("MRH")

MRH.default = function(x){
	if(nrow(x) == 1 | ncol(x) < 5){ 
		stop("Data set provided is not an MRH MCMC chain")
	}
	x = as.matrix(x)
	class(x) <- "MRH"
	x
}

is.MRH = function(x){	if(inherits(x, "MRH")){	TRUE } else {	FALSE	}	}
as.MRH = function(x){	UseMethod("MRH")	}

# print will eventually contain the call, etc (see lm results)
summary.MRH = function(object, alpha.level = .05,...){
	
	if(is.null(colnames(object))){ if(names(object)[1] == 'summary'){
			object = MRH(read.table(paste(object$outfolder, '/MCMCchains.txt', sep = ''), header = TRUE))
	}}
	if(substr(colnames(object)[2], 3, 3) == ""){	output = callPHsummary(object[,-ncol(object)], alpha.level = alpha.level)
	} else {	output = callBAsummary(object[,-ncol(object)], alpha.level = alpha.level)	}
	return(output)
}

plot.MRH = function(x, maxStudyTime, main = "", xlab = 'Time', ylab = 'Hazard Rate', plot.type = 'h', interval = TRUE,
alpha.level = 0.05, smooth.graph = FALSE, smooth.df = NULL, combine.graphs = TRUE, log.ratio = TRUE, ...){
	
	censortime = NULL
	if(!missing(maxStudyTime)){	censortime = maxStudyTime	}
	if(is.MRH(x) == FALSE){ stop("Object is not an MRH object")	}
	if(is.null(colnames(x))){ if(names(x)[1] == 'summary'){
		censortime = x$maxStudyTime
		x = MRH(read.table(paste(x$outfolder, '/MCMCchains.txt', sep = ''), header = TRUE))
	}} 
	if(is.null(censortime) & 'h' == plot.type){ stop("Maximum study time needed for hazard rate calculation")	}

	if(substr(colnames(x)[2], 3, 3) == ""){
		callPHplot(x, censortime, main = main, xlab = xlab, ylab = ylab, plot.type = plot.type, 
				   smooth.graph = smooth.graph, smooth.df = smooth.df, conf.int = interval, alpha.level = alpha.level)
	} else {
		callBAplot(x, censortime, main = main, xlab = xlab, ylab = ylab, plot.type = plot.type, 
				   smooth.graph = smooth.graph, smooth.df = smooth.df, combine.graphs = combine.graphs, 
				   conf.int = interval, alpha.level = alpha.level, log.ratio = log.ratio)
	}
}
