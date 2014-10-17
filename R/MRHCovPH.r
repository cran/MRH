###########################################################################################
# This code is used to estimate the parameters of a given data set using the MRH
#	methodology.  One to many covariates can be included under the proportional 
#	hazards assumption.
###########################################################################################

MRHCovPH = function(Mval, Ti, delta, X, outfilename, censortime,
prune.indc, burnIn, maxIter, thin, RmpInit, a.init, lambda.init, beta.init, 
k, gamma.mp, GR, convGraphs, fix.burnIn, fix.thin, fix.max,
writenum, sysseconds, systemtime, checknum){
	
	###########################################################################
	#		If X matrix is available, standardize the continuous 
	#		variables and label the betas
	###########################################################################
	if(!is.null(X)){
		n = length(Ti)
		betaNames = colnames(X)
		numParams = ncol(X)
				
		# Find the bounds on the beta estimates, used for sampling
		betaInfo = summary(coxph(Surv(Ti, delta) ~ as.matrix(X)))$coeff
		if(length(which(betaInfo[,1] < -50 | betaInfo[,1] > 50)) > 0){
			stop(paste("Algorithm will not converge.  Covariate estimate for variable(s)", 
					   which(betaInfo[,1] < -50 | betaInfo[,1] > 50), "is out of range."))
		}
		betaLB = betaInfo[,1]-betaInfo[,3]*5*sqrt(n)
		betaUB = betaInfo[,1]+betaInfo[,3]*5*sqrt(n)
		betaLB[betaLB < -50] = -50
		betaUB[betaUB > 50] = 50
		betaInfo[,1][betaInfo[,1] < betaLB | betaInfo[,1] > betaUB] = 0
			
		# Check the bounds on the continuous covariates 
		standardizeinfo = prepX.contin(X, betaLB, betaUB)
		Xstdz = as.matrix(standardizeinfo$X)
		stdIndex = standardizeinfo$stdIndex
		Xmeans = rep(0, numParams)
		Xmeans[stdIndex] = standardizeinfo$Xmeans
		Xsds = rep(1, numParams)
		Xsds[stdIndex] = standardizeinfo$Xsds
		X = as.matrix(X)
		
		if(length(stdIndex) > 0){
			betas.LB.stdz = betaLB*Xsds
			betas.UB.stdz = betaUB*Xsds
			betas.LB.stdz[betas.LB.stdz < -50] = -50
			betas.UB.stdz[betas.UB.stdz > 50] = 50
		}
	} else {
		numParams = 0
		n = length(Ti)
		Xmeans = Xsds = stdIndex = NULL
		X = matrix(1, ncol = 1, nrow = n)
		beta = matrix(0, ncol = 1, nrow = 1)
	}
		
	#########################################################
	# Create all variables that only need to be set up once 
	# during the estimation.
	#########################################################
	# mat01 is the PI matrix on p328 of the SINICA paper
	mat01 = calc_mat01(M = Mval)
	# TiBR holds the bin number for each persons time, as well
	# as the ratio of time each person spent in the final bin.	
	TiBR = getbin.and.ratio(TJ = censortime, Ti = Ti, M = Mval)

	# Calculate temp, which holds the indicators for both the Rmp and 1-Rmp values.
	# These indicators are needed to calculate the Rmp posteriors.
	temp = calcXIndics(mat01 = mat01, bin_ratio = TiBR, M = Mval)
	# RmpInd holds the indicator for the Rmp values, with the columns as follows: 
	# {R10, R20, R21, R30, etc}. The index goes from m = 1 to M, and p = 0 to 2^(m-1)-1.
	RmpInd = matrix(temp[,1:(2^Mval-1)], ncol = 2^Mval-1)
	# one_RmpInd holds the 1-Rmp indicators, with the columns as follows:
	# {1-R10, 1-R20, 1-R21, etc}. The index goes from m = 1 to M, and p = 0 to 2^(m-1)-1.
	one_RmpInd = matrix(temp[,-c(1:(2^Mval-1))], ncol = 2^Mval-1)
	
	# The inBin and failBin matrices are needed for calcluating H and F quickly.  
	whereFail = TiBR$bin+1
	whereFail[whereFail > 2^Mval] = 2^Mval
	# inBin is a matrix that holds a '1' for every bin the subject 
	# survived through, and the ratio in the bin they failed in, with
	# zeroes in all columns beyond the failure bin.
	# failBin is a matrix that holds a '1' in the column 
	# corresponding to the bin that the subject failed in, and '0' otherwise.
	inBin = failBin = matrix(0, ncol = 2^Mval, nrow = n)
	inBin[,1] = TiBR$ratio
	failBin[whereFail == 1, 1] = 1
	for(j in 2:2^Mval){
		inBin[whereFail == j, 1:(j-1)] = 1
		inBin[whereFail == j, j] = TiBR$ratio[whereFail == j]
		failBin[whereFail == j, j] = 1
	}
	
	# mvec will hold the value of m for m = 1,...,M, with 
	# each number repeated 2^(m-1) times.  This is used in 
	# the calculations of the Rmp posteriors.
	mvec = 1
	if(Mval > 1){ for(mctr in 2:Mval){mvec = c(mvec, rep(mctr, 2^(mctr-1)))}}
		
	# Make the labels for the Rmps so that the output is easier to read	
	RmpRowNums = NULL
	for(mctr in 1:Mval){ for(pctr in 0:(2^(mctr-1)-1)){ 
		RmpRowNums = c(RmpRowNums, paste(mctr, '.', pctr, sep = '')) 
	}}
	
	#########################################################
	# Initialize parameters
	#########################################################
	Rmp = rep(0.5, 2^Mval-1)
	a = a.init
	lambda = -log(mean(delta))/a
	if(!is.null(lambda.init)){	lambda = lambda.init	}
	
	# If user has specified their own Rmp initialization values, then use them
	if(!is.null(RmpInit)){ Rmp = RmpInit	}
	if(numParams > 0){
		if(is.null(beta.init)){	beta = betaInfo[,1]
		} else {	beta = beta.init	}
	} 
	# If the user is not using pruned bins, set Rmp sampling index vector so that 
	# all Rmps are sampled.  If the user would like to use pruned bins, set the 
	# Rmp sampling index vector so that those Rmps are not sampled or even checked.
	RmpNoPruneIndx = 1:(2^Mval-1)
	if(!is.null(prune.indc)){  RmpNoPruneIndx = RmpNoPruneIndx[-which(prune.indc == 1)]	}
	# If user is going to use the Gelman-Rubin test statistic, then jiggle the
	# initial values to cover the range for each parameter.
	if(GR == TRUE){
		Rmp[RmpNoPruneIndx] = runif(length(RmpNoPruneIndx), 0, 1)
		a = a + runif(1, 0, 10)
		lambda = lambda + runif(1, 0, 1)
		if(numParams > 0){	for(i in 1:numParams){	beta[i] = runif(1, betaLB[i], betaUB[i])	}}
	}
	# Record all initialized values
	tempRmp = as.data.frame(matrix(Rmp, nrow = 1))
	names(tempRmp) = paste('Rmp',RmpRowNums, sep = '')
	if(numParams > 0){
		tempbetas = as.data.frame(matrix(beta, nrow = 1))
		names(tempbetas) = betaNames	
	} else {	tempbetas = NULL	}
	initialValues = list(RmpInit = tempRmp, betaInit = tempbetas, aInit = a, lambdaInit = lambda)
	
	#########################################################
	# Start output file for parameter estimates
	#########################################################
	# Write out ds, betas, H00, Rmps
	if(numParams > 0){
		if(length(betaNames) != numParams){	betaNames = 1:numParams	}
		write.table(matrix(c('iteration', paste('d',1:2^Mval, sep = ''), paste('beta', betaNames, sep = '.'), 
				'H00', paste('Rmp', RmpRowNums, sep = '')), nrow = 1), 
				paste(outfilename, '/MCMCchains.txt', sep = ''), row.names = FALSE, col.names = FALSE)
		round.values = rep(16, length(c(paste('d',1:2^Mval, sep = ''), paste('beta', betaNames, sep = '.'), 
										'H00', paste('Rmp', RmpRowNums, sep = ''))))
	} else {
		write.table(matrix(c('iteration', paste('d',1:2^Mval, sep = ''),  
					'H00', paste('Rmp', RmpRowNums, sep = '')), nrow = 1), 
					paste(outfilename, '/MCMCchains.txt', sep = ''), row.names = FALSE, col.names = FALSE)
		round.values = rep(16, length(c(paste('d',1:2^Mval, sep = ''), 'H00', paste('Rmp', RmpRowNums, sep = ''))))
	}
	#########################################################
	# Run Gibbs sampler
	#########################################################
	convergence = best.thin.found = FALSE
	if(fix.thin == TRUE){ best.thin.found = TRUE }
	outdata = NULL
	loopctr = 1

	while((loopctr <= maxIter & convergence == FALSE & fix.max == FALSE) | 
		  (loopctr <= maxIter & fix.max == TRUE)){
		
		######## Draw H00 ########
		# Calculate F, which is needed in the posterior of H00
		F = calc_F(mat01 = mat01, Rmp = Rmp, M = Mval, inBinMat = inBin)
		# Draw H00
		H00 = rgamma(1, shape = a+sum(delta), scale = (1/lambda+sum(exp(X%*%beta)*F))^-1)
	
		######## Draw Rmp values ########
		# Calculate H, which is needed in the posterior of the Rmps
		H = calc_H(mat01 = mat01, H00 = H00, Rmp = Rmp, M = Mval, inBinMat = inBin)
		# Draw Rmps
		for(Rmpctr in RmpNoPruneIndx){
			Rmp[Rmpctr] = .Call("arms", c(.001, .999), f = function(x) logRmpPost_PHcovs(x, RmpFull = Rmp, H00val = H00, kval = k, 
																						 aval = a, 
																						 gamma.mpval = gamma.mp[Rmpctr],
																						 betavec = beta, X.mat = X, 
																						 mval = mvec[Rmpctr], 
																						 RmpIndic = RmpInd[,Rmpctr], 
																						 one_RmpIndic = one_RmpInd[,Rmpctr],
																						 deltavals = delta, Mvalue = Mval, 
																						 inBinMat = inBin, mat01 = mat01, 
																						 formula = Rmpctr), Rmp[Rmpctr], 
								as.integer(1), new.env(), PACKAGE = "MRH")
		}
		
		##### Draw beta ########
		# H needs to be recalculated first
		if(numParams > 0){
			if(length(stdIndex) > 0){
				betas.stdz = beta*Xsds
				H00.stdz = H00*exp(sum(betas.stdz/Xsds*Xmeans))
				H = calc_H(mat01 = mat01, H00 = H00.stdz, Rmp = Rmp, M = Mval, inBinMat = inBin)
				for(betaCtr in 1:numParams){
					H00.stdz = H00*exp(sum(betas.stdz/Xsds*Xmeans))
					H = calc_H(mat01 = mat01, H00 = H00.stdz, Rmp = Rmp, M = Mval, inBinMat = inBin)
					betas.stdz[betaCtr] = .Call("arms", c(betas.LB.stdz[betaCtr], betas.UB.stdz[betaCtr]), 
													   f = function(x) logbetaPost_PH(x, betaFull = betas.stdz, 
																					  deltavals = delta, Xmatrix = Xstdz, 
																					  Hvals = H, whichBeta = betaCtr, mu.beta = 0, 
																					  sigma.beta = 10), betas.stdz[betaCtr], 
													   as.integer(1), new.env(), PACKAGE = "MRH")
				}
				beta = betas.stdz/Xsds
			} else {
				H = calc_H(mat01 = mat01, H00 = H00, Rmp = Rmp, M = Mval, inBinMat = inBin)
				for(betaCtr in 1:numParams){
					beta[betaCtr] = .Call("arms", c(betaLB[betaCtr], betaUB[betaCtr]), 
										  f = function(x) logbetaPost_PH(x, betaFull = beta, deltavals = delta, 
																		 Xmatrix = X, Hvals = H, whichBeta = betaCtr, mu.beta = 0, 
																		 sigma.beta = 10), beta[betaCtr], 
										  as.integer(1), new.env(), PACKAGE = "MRH")
				}
			}
		}
		##### Draw lambda ########
		lambda = .Call("arms", c(.0001, 10), f = function(x) logLambdaPost(x, mu.lval = 100, H00val = H00, aval = a),
					   lambda, as.integer(1), new.env(), PACKAGE = "MRH")
		##### Draw a #########
		a = sample.a(1:50, lambdaval = lambda, mu.aval = 4, H00val = H00, kval = k, 
							 mvec = mvec, Mval = Mval, Rmpvec = Rmp, gamma.mpvec = gamma.mp)

		#######################################################################################
		#					Gather and provide information for user 
		#######################################################################################
		# If it is every thin_th iteration, save to the data set
		if(loopctr == 100){	
			currenttime = unclass(Sys.time())
			estruntime = round((currenttime-unclass(systemtime)-.05)*min(c(maxIter, 100000))/100)
			units = 'seconds'
			if(estruntime > 60){ 
				units = 'minutes'
				estruntime = round(estruntime/60)
			}
			if(estruntime > 60){ 
				units = 'hours'
				estruntime = ceiling(estruntime/60)
			}
			writeLines(paste('Estimated total run time is', estruntime, units, '\n'))
			writeLines('To shorten the run time, re-run with fewer iterations or a smaller number of bins. \n')
		}
		if((loopctr %% thin) == 0 & loopctr >= burnIn){ 
			d.iter = calc_dfast(Mval, Rmp, H00, mat01)
			if(numParams > 0){
				outdata = rbind(outdata, matrix(c(loopctr, d.iter, beta, H00, Rmp), nrow = 1))
			} else {
				outdata = rbind(outdata, matrix(c(loopctr, d.iter, H00, Rmp), nrow = 1))
			}
		}

		# Every 10,000 or writenum iterations, write out the data set and empty the working one.
		if((loopctr %% writenum) == 0 & loopctr >= burnIn){
			if(nrow(outdata) > 1){
				outdata = cbind(outdata[,1], sapply(2:ncol(outdata), function(x) round(outdata[,x], round.values[x-1])))
			} else {
				outdata = matrix(c(outdata[,1], sapply(2:ncol(outdata), function(x) round(outdata[,x], round.values[x-1]))), nrow - 1)
			}
			write.table(outdata, paste(outfilename, '/MCMCchains.txt', sep = ''), 
						col.names = FALSE, row.names = FALSE, append = TRUE)
			if(loopctr != maxIter){	outdata = NULL	}
		}
		# Every 100,000 or checknum iterations, check for convergence with the full data set (minus the burned values)
		if((loopctr %% checknum) == 0 & loopctr >= burnIn){
			testData = read.table(paste(outfilename, '/MCMCchains.txt', sep = ''), header = TRUE)
			
			if(loopctr == 100000){	round.values = FindRoundVals(testData)	}
				
			# Test autocorrelation and find best thinning value 
			while(best.thin.found == FALSE){
				newthin = thinAmt(thin, mcmc(testData[,c((2^Mval+2):(2^Mval+2+numParams),
										((2^Mval+3+numParams):(2^(Mval+1)+numParams+1))[RmpNoPruneIndx])], 
											 thin = thin, start = burnIn))
				if(newthin == thin){ 
					best.thin.found = TRUE
				} else {	
					testData = testData[seq(1, nrow(testData), by = newthin/thin),]
					thin = newthin	
				}
			}
			keepnames = colnames(testData)
			testData = cbind(testData[,1], sapply(2:ncol(testData), function(x) round(testData[,x], round.values[x-1])))
			colnames(testData) = keepnames
			write.table(testData, paste(outfilename, '/MCMCchains.txt', sep = ''), row.names = FALSE)
			
			# Test convergence
			enough.data.left = TRUE
			if(nrow(testData) < 500){ enough.data.left = FALSE	}
			newburnIn = burnIn
			while(convergence == FALSE & enough.data.left == TRUE){
				convergence = convergeFxn(mcmc(testData[,c((2^Mval+2):(2^Mval+2+numParams),
								((2^Mval+3+numParams):(2^(Mval+1)+numParams+1))[RmpNoPruneIndx])], 
											   thin = thin, start = burnIn))
				if(fix.burnIn == FALSE & convergence == FALSE & nrow(testData) > 30000/thin){
						testData = testData[-c(1:(20000/thin)),]
						newburnIn = newburnIn + 20000
				} else { enough.data.left = FALSE	}
			}
			if(convergence == TRUE & fix.burnIn == TRUE){
				burnIn = newburnIn
				keepnames = colnames(testData)
				testData = cbind(testData[,1], sapply(2:ncol(testData), function(x) round(testData[,x], round.values[x-1])))
				colnames(testData) = keepnames
				write.table(testData, paste(outfilename, '/MCMCchains.txt', sep = ''), row.names = FALSE)
			}
		}
		if((loopctr %% 5000) == 0){	print(paste(loopctr, 'iterations completed...'))	}
		
		# Increment loopctr
		loopctr = loopctr+1
	}

	print("Estimation routine finished, preparing results....", quote = FALSE)
	
	#############################################################################################################
	# Make convergence graphs, hazard rate graphs (raw and smoothed), calculate summaries, calculate 
	# AIC/BIC/DIC, add log-likelihood to MCMCchains file, and return to user.
	#############################################################################################################
	finaldata = read.table(paste(outfilename, '/MCMCchains.txt', sep = ''), header = TRUE)
	# Convergence graphs #
	if(convGraphs == TRUE){	
		convergenceGraphs(finaldata[,c((2^Mval+2):(2^Mval+2+numParams),
									   ((2^Mval+3+numParams):(2^(Mval+1)+numParams+1))[RmpNoPruneIndx])], 
						  paste(outfilename, '/convergenceGraphs', sep = ''), burnIn, thin)
	}
	# Summaries of MCMCchains	
	outputdata = AnalyzeComplete(Mval, censortime, finaldata[,-1], paste(outfilename, '/', sep = ''), betaNames)
	# Calculation of AIC, BIC, DIC and add likelihood to chains
	nrow.data = nrow(finaldata)
	binwidth = censortime/2^Mval
	Rmpvec = as.data.frame(matrix(NA, ncol = sum(2^(1:Mval)), nrow = nrow.data))
	Rmpvec[,1:(sum(2^(1:Mval))/2)*2-1] = finaldata[,1:(2^Mval-1)+2^Mval+1+numParams+1]
	Rmpvec[,1:(sum(2^(1:Mval))/2)*2] = 1-finaldata[,1:(2^Mval-1)+2^Mval+1+numParams+1]
	logds = matrix(rep(log(finaldata$H00), each = 2^Mval), ncol = 2^Mval, nrow = nrow.data, byrow = TRUE) + 
							t(mat01%*%t(log(Rmpvec))) - log(binwidth)
	if(numParams > 0){	Xbeta = rowSums(as.matrix(X[rep(1:n, nrow.data),])*(as.matrix(finaldata[,1:numParams + 2^Mval+1])[rep(1:nrow.data, each = n),]))
	} else {	Xbeta = rep(0, n*nrow.data)	}
	totalone = rep(delta, nrow.data)*(rowSums(failBin[rep(1:n, nrow.data),]*(logds[rep(1:nrow.data, each = n),]))+Xbeta) - 
				rowSums(inBin[rep(1:n, nrow.data),]*(finaldata[,1:2^Mval+1][rep(1:nrow.data, each = n),]))*exp(Xbeta)
	logliks = rowSums(matrix(totalone, ncol = n, byrow = TRUE))
	finaldata = as.data.frame(cbind(finaldata, logliks))
	names(finaldata) = c(names(finaldata)[-ncol(finaldata)], 'loglike')
	write.table(finaldata, paste(outfilename, '/MCMCchains.txt', sep = ''), row.names = FALSE)
	# Number of estimated parameters equals (a, lambda) + H00 + number of estimated Rmps + number of covariates
	numberofEstParams = 2 + 1 + numParams + length(RmpNoPruneIndx)
	AIC = 2*(numberofEstParams) - 2*min(logliks)
	BIC = -2*min(logliks) + (numberofEstParams)*log(n)
	DIC = .5*var(-2*logliks) + mean(-2*logliks)
	
	hazardRate = outputdata$d/(censortime/2^Mval)
	colnames(hazardRate) = paste('q', c('.5', '.025', '.975'), sep = '')
	rownames(hazardRate) = paste('h.Bin', 1:2^Mval, sep = '')
	
	if(GR == FALSE){ gr.used = 'Option not used' } else { gr.used = 'Option Used'	}
	outputdata = c(outputdata, list(hazardRate = hazardRate, AIC = AIC, BIC = BIC, DIC = DIC, 
									burnIn = burnIn, thin = thin, TotalIters = loopctr-1, 
									convergence = convergence, gelman.rubin.used = gr.used, fix.thin = fix.thin,
									fix.burnIn = fix.burnIn, fix.max = fix.max, initialValues = initialValues,
									runtime = paste(round((unclass(Sys.time()) - unclass(systemtime))/(60*60), 2), 'hours'), 
									outfolder = outfilename, maxStudyTime = censortime))
	class(outputdata) = "MRH"
	return(outputdata)
	
}
