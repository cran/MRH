###########################################################################################
# Contains the functions needed to: convert to an MRH object, summary, and plot for PH model.
###########################################################################################

callPHsummary = function(object, alpha.level){

	M = log(length(which(substr(colnames(object), 1, 1) == 'd')))/log(2)
	object = object[,-1]

	substr.4 = substr(colnames(object), 1, 4)
	substr.3 = substr(colnames(object), 1, 3)
	substr.1 = substr(colnames(object), 1, 1)
	any.betas = any.ks = FALSE
	if(length(which(substr.4 == 'beta')) > 0){	
		any.betas = TRUE	
		beta.loc = which(substr.4 == 'beta')
	}
	if(length(which(substr.1 == 'k')) > 0){	
		any.ks = TRUE	
		k.loc = ncol(object)
	}
	H.loc = which(substr.1 == 'H')
	Rmp.loc = which(substr.3 == 'Rmp')
	
	title.lb = alpha.level/2*1000
	if(alpha.level/2*1000 < 100){	title.lb = paste('0', title.lb, sep = '')	}
	if(alpha.level/2*1000 < 10){	title.lb = paste('0', title.lb, sep = '')	}
	title.ub = (1-alpha.level/2)*1000
	
	################ Calculate the outputs for each parameter ################
	### Calculate d
	dInts = t(sapply(1:2^M, function(x) quantile(object[,x], probs = c(.5, alpha.level/2, 1-alpha.level/2))))
	doutput = as.data.frame(dInts)
	names(doutput) = c('dEst', paste('dq', title.lb, sep = '.'), paste('dq', title.ub, sep = '.'))
	row.names(doutput) = paste('d', 1:2^M, sep = '')

	### Calculate beta if there are any ###
	if(any.betas == TRUE){
		betas = t(sapply(beta.loc, function(x) 
						 quantile(object[,x], probs = c(.5, alpha.level/2, 1-alpha.level/2))))
		betaoutput = as.data.frame(betas)
		betanames = colnames(object)[beta.loc]
		names(betaoutput) = c('betaEst', paste('betaq', title.lb, sep = '.'), paste('betaq', title.ub, sep = '.'))
		row.names(betaoutput) = betanames
	}
	if(any.ks == TRUE){
		k.output = as.data.frame(t(quantile(object[,k.loc], probs = c(.5, .025, .975))))
		names(k.output) = c('kEst', paste('kq', title.lb, sep = '.'), paste('kq', title.ub, sep = '.'))
		row.names(k.output) = 'k'
	}

	#### Calculate the Rmps 
	RmpInts = as.data.frame(t(sapply(Rmp.loc, 
							function(x) quantile(object[,x], probs = c(.5, alpha.level/2, 1-alpha.level/2)))))
	names(RmpInts) = c('RmpEst', paste('Rmpq', title.lb, sep = '.'), paste('Rmpq', title.ub, sep = '.'))
	RmpRowNums = NULL
	for(mctr in 1:M){ for(pctr in 0:(2^(mctr-1)-1)){ 
		RmpRowNums = c(RmpRowNums, paste(mctr, '_', pctr, sep = '')) 
	}}
	row.names(RmpInts) = paste('Rmp', RmpRowNums, sep = '')

	#### Calculate the cumulative hazard H
	Houtput = as.data.frame(t(quantile(object[,H.loc], probs = c(.5, alpha.level/2, 1-alpha.level/2))))
	names(Houtput) = c('HEst', paste('Hq', title.lb, sep = '.'), paste('Hq', title.ub, sep = '.'))
	row.names(Houtput) = c('H')

	if(any.betas == TRUE){	
		output = list(doutput, betaoutput, Houtput, RmpInts)
		names(output) = c('d', 'beta', 'H', 'Rmp')
	} else { 
		output = list(doutput, Houtput, RmpInts)
		names(output) = c('d', 'H', 'Rmp')
	}
	if(any.ks == TRUE){
		output = c(output, list(k = k.output))
	}

	return(output)
}

callPHplot = function(x, censortime, main, xlab, ylab, plot.type, 
smooth.graph, smooth.df, combine.graphs, conf.int, alpha.level){
	
	# Calculate M 
	M = log(length(which(substr(colnames(x), 1, 1) == 'd')))/log(2)
	# Cannot smooth with 2 bins, so set smooth = F if user has selected this option 
	# (the graph will still appear smooth to the user)
	if(M == 1 & smooth.graph == TRUE){	smooth.graph = FALSE	}

	# Get the parameter estimates #
	if(!missing(censortime)){
		results = t(CalcFunction(x, function.type = plot.type, alpha.level = alpha.level, maxStudyTime = censortime)[[1]])
	} else {	
		results = t(CalcFunction(x, function.type = plot.type, alpha.level = alpha.level)[[1]])
		censortime = 1
	}
	results = results[,c(2,1,3)]
	
	if(plot.type == 'H'){	ylab = 'Cumulative hazard'
	} else if(plot.type == 'S'){	ylab = 'Survival Function'	}
	
	###### Plot #####
	if(smooth.graph == FALSE){
		if(plot.type == 'h'){	xvals = pbfxn(2^M, censortime/2^M, results[,1])
		} else {
			xvals = pbfxn(2^M, censortime/2^M, results[-1,1])
			if(plot.type == 'H'){	xvals = rbind(c(0, 0), xvals)
			} else { xvals = rbind(c(0, 1), xvals)	}
		}
		ylimit = range(results)
		if(conf.int == FALSE){	ylimit = range(results)	}
		plot(xvals$x, xvals$y, main = main, xlab = xlab, ylab = ylab, type = 'l', 
			 ylim = range(results), lwd = 3)
		if(conf.int == TRUE){
			if(plot.type == 'h'){
				points(xvals$x, pbfxn(2^M, censortime/2^M, results[,2])$y, type = 'l', lty = 2)
				points(xvals$x, pbfxn(2^M, censortime/2^M, results[,3])$y, type = 'l', lty = 2)
			} else if(plot.type == 'H'){
				points(xvals$x, c(0, pbfxn(2^M, censortime/2^M, results[-1,2])$y), type = 'l', lty = 2)
				points(xvals$x, c(0, pbfxn(2^M, censortime/2^M, results[-1,3])$y), type = 'l', lty = 2)
			} else {
				points(xvals$x, c(1, pbfxn(2^M, censortime/2^M, results[-1,2])$y), type = 'l', lty = 2)
				points(xvals$x, c(1, pbfxn(2^M, censortime/2^M, results[-1,3])$y), type = 'l', lty = 2)
			}				
		}
	} else {	
		# Calculate and plot the smooth ds over time 
		if(is.null(smooth.df)){	smooth.df = 2^M/2	}
		timePts = seq(censortime/2^M, censortime, by = censortime/2^M)-.5*censortime/2^M
		if(plot.type != 'h'){	timePts = c(0, timePts)	}
		smoothdInts = smooth.spline(results[,1]~timePts, df = smooth.df)$y
		smoothdInts = cbind(smoothdInts, smooth.spline(results[,2]~timePts, df = smooth.df)$y)
		smoothdInts = cbind(smoothdInts, smooth.spline(results[,3]~timePts, df = smooth.df)$y)
		ylimit = range(smoothdInts)
		if(conf.int == FALSE){	ylimit = range(smoothdInts[,1])	}
		plot(timePts, smoothdInts[,1], main = main, xlab = xlab, ylab = ylab, type = 'l', 
			 ylim = range(smoothdInts), lwd = 3)
		if(conf.int == TRUE){
			points(timePts, smoothdInts[,2], type = 'l', lty = 2)
			points(timePts, smoothdInts[,3], type = 'l', lty = 2)
		}
	}
}	
