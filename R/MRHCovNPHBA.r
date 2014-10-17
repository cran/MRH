###########################################################################################
# Created 6Sep13
###########################################################################################

MRHCovNPHBA = function(Mval, Ti, Xfixed, XNPH, delta, outfilename, 
prune.indc, burnIn, maxIter, thin, RmpInit, a.init, k, gamma.mp, 
lambda.init, beta.init, censortime, GR, convGraphs, fix.burnIn, fix.thin, fix.max, 
writenum, sysseconds, systemtime, checknum){

	###########################################################################
	#	Standardize the continuous variables and label the betas
	###########################################################################
	n = length(Ti)
			
	# Convert categorical variables to 0/1 
	namesHazGroups = names(table(XNPH))
	betaNPHnames_preBin = paste(colnames(XNPH), namesHazGroups[-1], sep = '.') 
			
	# Find the number of PH and NPH parameters, and create the beta names	
	if(!is.null(Xfixed)){	numPHParams = ncol(Xfixed)	} else {	numPHParams = 0	}
	numNPHParams = ncol(XNPH)
	numParams = numPHParams + numNPHParams
	if(is.na(betaNPHnames_preBin)[1]){ betaNPHnames_preBin = paste(rep('var', numNPHParams), 1:numNPHParams, sep = '') }
	betaNames = paste(rep(paste(betaNPHnames_preBin, '.bin', sep = ''), each = 2^Mval), 1:2^Mval, sep = '')
	if(numPHParams > 0){	betaNames = c(colnames(Xfixed), betaNames)		}
	# Find the bounds on the beta estimates, used for sampling
	binFail = getbin.and.ratio(TJ = censortime, Ti = Ti, M = Mval)$bin + 1
	binFail[binFail > 2^Mval] = 2^Mval
	binvar = factor(binFail, levels = 1:2^Mval)
	if(numPHParams > 0){
		betaInfo = summary(coxph(Surv(Ti, delta) ~ as.matrix(Xfixed)))$coeff
		if(length(which(betaInfo[,1] < -25 | betaInfo[,1] > 25)) > 0){
			stop(paste("Algorithm will not converge.  Covariate estimate for variable(s)", 
					   betaNames[which(betaInfo[,1] < -25 | betaInfo[,1] > 25)], "is (are) out of range."))
		}
		betaLB = betaInfo[,1]-betaInfo[,3]*3*sqrt(n)
		betaUB = betaInfo[,1]+betaInfo[,3]*3*sqrt(n)
		betaLB[betaLB < -25 | is.na(betaInfo[,4])] = -25
		betaLB[betaLB > 0] = -1
		betaUB[betaUB > 25 | is.na(betaInfo[,4])] = 25
		betaUB[betaUB < 0] = 1
		betaInfo[is.na(betaInfo[,4]),1] = 0
	}
	# Check the bounds on the continuous covariates 
	if(numPHParams > 0){
		stdIndex = NULL
		Xmeans = rep(0, numPHParams)
		Xsds = rep(1, numPHParams)
		standardizeinfo = prepX.contin(Xfixed, betaLB[1:numPHParams], betaUB[1:numPHParams])
		Xstdz = as.matrix(standardizeinfo$X)
		stdIndex = standardizeinfo$stdIndex
		Xmeans[stdIndex] = standardizeinfo$Xmeans
		Xsds[stdIndex] = standardizeinfo$Xsds
	}
			
	#########################################################
	# Create all variables that only need to be set up once 
	# during the estimation.
	#########################################################
	# mat01 is the PI matrix on p328 of the SINICA paper
	mat01 = calc_mat01(M = Mval)
	# Holds the sum of the columns of mat01 in pairs of two, and is used in the Rmpj posterior.
	newmat01 = cbind(1, mat01[,1:(ncol(mat01)/2-1)])
	# TiBR holds the bin number for each persons time, as well
	# as the ratio of time each person spent in the final bin.	
	TiBR = getbin.and.ratio(TJ = censortime, Ti = Ti, M = Mval)
			
	# The inBin and failBin matrices are needed for estimating betaj parameters.  
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
			
	# Calculate temp, which holds the indicators for both the Rmp and 1-Rmp values.
	# These indicators are needed to calculate the Rmp posteriors.
	temp = calcXIndics(mat01 = mat01, bin_ratio = TiBR, M = Mval)
	# RmpInd holds the indicator for the Rmp values, with the columns as follows: 
	# {R10, R20, R21, R30, etc}. The index goes from m = 1 to M, and p = 0 to 2^(m-1)-1.
	RmpInd = matrix(temp[,1:(2^Mval-1)], ncol = 2^Mval-1)
	# one_RmpInd holds the 1-Rmp indicators, with the columns as follows:
	# {1-R10, 1-R20, 1-R21, etc}. The index goes from m = 1 to M, and p = 0 to 2^(m-1)-1.
	one_RmpInd = matrix(temp[,-c(1:(2^Mval-1))], ncol = 2^Mval-1)
	# mvec will hold the value of m for m = 1,...,M, with 
	# each number repeated 2^(m-1) times.  This is used in 
	# the calculations of the Rmp posteriors.
	mvec = 1
	if(Mval > 1){for(mctr in 2:Mval){	mvec = c(mvec, rep(mctr, 2^(mctr-1)))	}}
	# ind and one_ind are also needed in the calculation of the Rmp posteriors.
	ind = matrix(mat01[,1:(sum(2^(1:Mval))/2)*2-1], ncol = 2^Mval-1)
	one_ind = matrix(mat01[,1:(sum(2^(1:Mval))/2)*2], ncol = 2^Mval-1)
			
	# IDtypes will hold the types of different hazards #
	IDtypes = names(table(XNPH))
	IDtypes = IDtypes[order(IDtypes)]
	numHazards = length(IDtypes)
	indices = list(which(XNPH == IDtypes[1]))
	for(hazCtr in 2:numHazards){
		indices = c(indices, list(which(XNPH == IDtypes[hazCtr])))
	}
		
	#########################################################
	# Initialize parameters
	#########################################################
	Rmp = matrix(0.5, nrow = 2^Mval-1, ncol = numHazards)
	a = a.init
	if(length(a.init) < numHazards){	a = rep(a.init, numHazards)	}
	if(!is.null(lambda.init)){	
		if(length(lambda.init) < numHazards){	lambda = rep(lambda.init, numHazards)	
		} else {	lambda = lambda.init	}
	} else {	lambda = -log(sapply(1:numHazards, function(x) mean(delta[indices[[x]]]))[])/a	}
	# If user has specified their own Rmp initialization values, then use them
	if(!is.null(RmpInit)){ Rmp = RmpInit	}
	# Initialize the beta estimates
	if(is.null(beta.init) & !is.null(Xfixed)){	betas = betaInfo[,1]
	} else if(!is.null(beta.init)){	betas = beta.init	
	} else {	
		betas = 0	
		Xfixed = matrix(1, ncol = 1, nrow = n)
	}
	# If the user is not using pruned bins, set Rmp sampling index vector so that 
	# all Rmps are sampled.  If the user would like to use pruned bins, set the 
	# Rmp sampling index vector so that those Rmps are not sampled or even checked.
	RmpNoPruneIndx = NULL
	numRmpsEstimated = 0
	for(rmpctr in 1:numHazards){
		numRmpsEstimated = numRmpsEstimated + 2^Mval-1
		if(!is.null(prune.indc)){  
			RmpNoPruneIndx = c(RmpNoPruneIndx, list((1:(2^Mval-1))[which(prune.indc[,rmpctr] == 0)]))	
			numRmpsEstimated = numRmpsEstimated - length(which(prune.indc[,rmpctr] == 1))
		} else {	RmpNoPruneIndx = c(RmpNoPruneIndx, list(1:(2^Mval-1)))	}
	}
	
	# If user is going to use the Gelman-Rubin test statistic, then jiggle the
	# initial values to cover the range for each parameter.
	if(GR == TRUE){
		for(rmpctr in 1:numHazards){
			Rmp[RmpNoPruneIndx[[rmpctr]],rmpctr] = runif(length(RmpNoPruneIndx[[rmpctr]]), 0, 1)
		}
		a = runif(numHazards, 5, 15)
		lambda = runif(numHazards, 0, .1)
		for(i in 1:length(betas)){	betas[i] = betas[i] + rnorm(1, mean = 0, sd = .25)	}
	}
			
	#########################################################
	# Start output file for parameter estimates
	#########################################################
	# Make the labels for the Rmps so that the output is easier to read	
	RmpRowNums = NULL
	for(mctr in 1:Mval){ for(pctr in 0:(2^(mctr-1)-1)){ 
		RmpRowNums = c(RmpRowNums, paste(mctr, '.', pctr, sep = '')) 
	}}
	RmpRowNums = paste(RmpRowNums, rep(namesHazGroups, each = 2^Mval-1), sep = '_')
	dnames = paste(paste('d', 1:2^Mval, sep = ''), rep(namesHazGroups, each = 2^Mval), sep = '_')
	write.table(matrix(c('iteration', dnames, paste('beta.', betaNames, sep = ''),
						 paste('H00', namesHazGroups, sep = '_'), paste('Rmp', RmpRowNums, sep = '')), nrow = 1), 
						 paste(outfilename, '/MCMCchains.txt', sep = ''), row.names = FALSE, col.names = FALSE)
	round.values = rep(16, length(c(dnames, paste('beta.', betaNames, sep = ''),
							 paste('H00', namesHazGroups, sep = '_'), paste('Rmp', RmpRowNums, sep = ''))))

	# Record all initialized values
	tempRmpj = as.data.frame(matrix(Rmp, nrow = 1))
	names(tempRmpj) = paste('Rmp',RmpRowNums, sep = '')
	if(!is.null(Xfixed)){	
		tempbetas = as.data.frame(matrix(betas, nrow = 1))	
		names(tempbetas) = betaNames[1:numPHParams]
	} else {	tempbetas = NULL	}
	initialValues = list(RmpInit = tempRmpj, betaInit = tempbetas, aInit = a, lambdaInit = lambda)	
			
	# analyzeIndex holds the columns that need to be analyzed in the burn-in, thinning, and 
	# convergence diagnosts/graphs.  This is just easier than writing out the indices every time those
	# functions are used.
	RmpAnalyzeIndex = 1:((2^Mval-1)*numHazards)
	if(!is.null(prune.indc)){	
		RmpAnalyzeIndex = rep(NA, (2^Mval-1)*numHazards)
		for(rmpctr in 1:numHazards){		
			RmpAnalyzeIndex[1:(2^Mval-1) + (rmpctr-1)*(2^Mval-1)][which(prune.indc[,rmpctr] == 0)] = 
				which(prune.indc[,rmpctr] == 0)+(rmpctr-1)*(2^Mval-1)
		}
	}
	RmpAnalyzeIndex = RmpAnalyzeIndex[!is.na(RmpAnalyzeIndex)]
	analyzeIndex = c(1:numHazards + 2^Mval*(2*numHazards-1) + numPHParams + 1, 
					 RmpAnalyzeIndex + 2^Mval*(2*numHazards-1) + numPHParams + numHazards + 1)
	if(numPHParams > 0){	analyzeIndex = c(1:numPHParams + 2^Mval*numHazards + 1, analyzeIndex)	}
	
	#########################################################
	# Run Gibbs sampler
	#########################################################
	outdata = NULL
	loopctr = 1
	convergence = best.thin.found = FALSE
	F = H = rep(NA, length = length(Ti))
	if(fix.thin == TRUE){ best.thin.found = TRUE }
	while((loopctr <= maxIter & convergence == FALSE & fix.max == FALSE) | 
		  (loopctr <= maxIter & fix.max == TRUE)){

		######## Draw H00 ########
		# Calculate F for each hazard group, which is needed in the posterior of each H00
		for(hazCtr in 1:numHazards){
			F[indices[[hazCtr]]] = calc_F(mat01 = mat01, Rmp = Rmp[,hazCtr], M = Mval, 
										  inBinMat = inBin[indices[[hazCtr]],])
		}			
		# Draw H00
		H00 = NULL
		for(hazCtr in 1:numHazards){
			H00 = c(H00, rgamma(1, shape = a[hazCtr]+sum(delta[indices[[hazCtr]]]), 
					scale = (1/lambda[hazCtr]+sum(exp(as.matrix(Xfixed[indices[[hazCtr]],])%*%betas)*
												  F[indices[[hazCtr]]]))^-1))
		}
		######## Draw Rmp values ########
		# Calculate H, which is needed in the posterior of the Rmps
		for(hazCtr in 1:numHazards){
			H[indices[[hazCtr]]] = calc_H(mat01 = mat01, H00 = H00[hazCtr], Rmp = Rmp[,hazCtr], M = Mval, 
							inBinMat = inBin[indices[[hazCtr]],])
		}
		# Draw Rmps
		for(hazCtr in 1:numHazards){
			for(rmpctr in RmpNoPruneIndx[[hazCtr]]){
				Rmp[rmpctr, hazCtr]= .Call("arms", c(.001, .999), f = function(x) logRmpPost_nonPHBA(x, 
						RmpFull = Rmp[,hazCtr], H00val = H00[hazCtr], kval = k, aval = a[hazCtr], 
						gamma.mpval = gamma.mp[rmpctr,hazCtr],
						betavec = as.matrix(betas), X.mat = as.matrix(Xfixed[indices[[hazCtr]],]), 
						mval = mvec[rmpctr], RmpIndic = RmpInd[indices[[hazCtr]],rmpctr], 
						one_RmpIndic = one_RmpInd[indices[[hazCtr]],rmpctr],
						deltavals = delta[indices[[hazCtr]]], Mvalue = Mval, 
						inBinMat = inBin[indices[[hazCtr]],], mat01 = mat01, 
						formula = rmpctr), Rmp[rmpctr,hazCtr], as.integer(1), new.env(), PACKAGE = "MRH")
			}
		}
		
		##### Draw beta values ########
		# H needs to be recalculated first
		for(hazCtr in 1:numHazards){
			H[indices[[hazCtr]]] = calc_H(mat01 = mat01, H00 = H00[hazCtr], 
										  Rmp = Rmp[,hazCtr], M = Mval, inBinMat = inBin[indices[[hazCtr]],])
		}
		if(numPHParams > 0){
			if(length(stdIndex) > 0){
				betas.stdz = betas*Xsds
				H00.stdz = H00*exp(sum(betas.stdz/Xsds*Xmeans))
				for(hazCtr in 1:numHazards){
					H[indices[[hazCtr]]] = calc_H(mat01 = mat01, H00 = H00.stdz[hazCtr], Rmp = Rmp[,hazCtr], M = Mval, 
												  inBinMat = inBin[indices[[hazCtr]],])
				}
				for(betaCtr in 1:numPHParams){
					betas.stdz = betas*Xsds
					H00.stdz = H00*exp(sum(betas.stdz/Xsds*Xmeans))
					for(hazCtr in 1:numHazards){
						H[indices[[hazCtr]]] = calc_H(mat01 = mat01, H00 = H00.stdz[hazCtr], Rmp = Rmp[,hazCtr], M = Mval, 
													  inBinMat = inBin[indices[[hazCtr]],])
					}
					betas.stdz[betaCtr] = .Call("arms", c(betaLB[betaCtr], betaUB[betaCtr]), f = function(x) 
												logbetaPost_NPHBA(x, betaFull = betas.stdz, deltavals = delta, Xmatrix = Xstdz, Hvals = H, 
																  whichBeta = betaCtr, mu.beta = 0, sigma.beta = 10), betas.stdz[betaCtr], 
												as.integer(1), new.env(), PACKAGE = "MRH")
				}
				betas = betas.stdz/Xsds
			} else {
				for(betaCtr in 1:ncol(Xfixed)){
					betas[betaCtr] = .Call("arms", c(betaLB[betaCtr], betaUB[betaCtr]), f = function(x) 
										   logbetaPost_NPHBA(x, betaFull = betas, deltavals = delta, Xmatrix = Xfixed, Hvals = H, 
														 whichBeta = betaCtr, mu.beta = 0, sigma.beta = 10), betas[betaCtr], 
										   as.integer(1), new.env(), PACKAGE = "MRH")
				}
			}
		}
		##### Draw lambda ########
		for(hazCtr in 1:numHazards){
			lambda[hazCtr] = .Call("arms", c(.0001, 10), f = function(x) logLambdaPost(x, 
					mu.lval = 100, H00val = H00[hazCtr], aval = a[hazCtr]), lambda[hazCtr], as.integer(1), new.env(), PACKAGE = "MRH")
		}
		
		##### Draw a #########
		for(hazCtr in 1:numHazards){
			a[hazCtr] = sample.a(1:50, lambdaval = lambda[hazCtr], mu.aval = 4, 
								 H00val = H00[hazCtr], kval = k, mvec = mvec, Mval = Mval, 
								 Rmpvec = Rmp[,hazCtr], gamma.mpvec = gamma.mp[,hazCtr])
		}
		
		#######################################################################################
		#					Gather and provide information for user 
		#######################################################################################
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
		# If it is every thin_th iteration, save to the data set
		if((loopctr %% thin) == 0 & loopctr >= burnIn){ 
			d.iter = betaholder = NULL
			for(hazCtr in 1:numHazards){
				d.iter = c(d.iter, calc_dfast(Mval, Rmp[,hazCtr], H00[hazCtr], mat01))
				if(hazCtr > 1){
					betaholder = c(betaholder, log(d.iter[2^Mval*(hazCtr-1)+1:2^Mval]/d.iter[1:2^Mval]))
				}
			}
			if(numPHParams > 0){
				outdata = rbind(outdata, 
							matrix(c(loopctr, d.iter, betas, betaholder,
								H00, matrix(Rmp, nrow = 1, byrow = TRUE)), nrow = 1))
			} else {
				outdata = rbind(outdata, 
								matrix(c(loopctr, d.iter, betaholder,
										 H00, matrix(Rmp, nrow = 1, byrow = TRUE)), nrow = 1))
			}	
		} 
		
		# Every 10,000 or writenum iterations, write out the data set and empty the working one.
		if((loopctr %% writenum) == 0 & loopctr >= burnIn){
			if(nrow(outdata) > 1){
				outdata = cbind(outdata[,1], sapply(2:ncol(outdata), function(x) round(outdata[,x], round.values[x-1])))
			} else {
				outdata = matrix(c(outdata[,1], sapply(2:ncol(outdata), function(x) round(outdata[,x], round.values[x-1]))), nrow = 1)
			}
			write.table(outdata, paste(outfilename, '/MCMCchains.txt', sep = ''), 
						col.names = FALSE, row.names = FALSE, append = TRUE)
			if(loopctr != maxIter){	outdata = NULL	}
		}
		# Every 100,000 or checknum iterations, check for convergence with the full data set (minus the burned values)
		if((loopctr %% checknum) == 0 & loopctr >= burnIn){
			testData = read.table(paste(outfilename, '/MCMCchains.txt', sep = ''), header = TRUE)
			
			# Find the right values for rounding
			if(loopctr == 100000){	round.values = FindRoundVals(testData)	}
			# Test autocorrelation and find best thinning value		
			while(best.thin.found == FALSE){
				newthin = thinAmt(thin, mcmc(testData[,analyzeIndex], thin = thin, start = burnIn))
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
			if(nrow(testData) < 500){	enough.data.left = FALSE	}
			newburnIn = burnIn
			while(convergence == FALSE & enough.data.left == TRUE){
				convergence = convergeFxn(mcmc(testData[,analyzeIndex], thin = thin, start = burnIn))
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
	# Convergence graphs
	if(convGraphs == TRUE){		
		convergenceGraphs(finaldata[,analyzeIndex], paste(outfilename, '/convergenceGraphs', sep = ''), burnIn, thin)
	}
	# Hazard rate graphs and summaries of parameter chains
	for(hazCtr in 1:numHazards){
		if(hazCtr == 1){
			outputdata = AnalyzeComplete(Mval, censortime, 
				finaldata[,c(1:2^Mval+2^Mval*(hazCtr-1) + 1,
					(2^Mval*numHazards+2):(2^Mval*(2*numHazards-1)+1+numPHParams),
					2^Mval*(2*numHazards-1)+numPHParams+1+hazCtr, 
					(2^Mval*(2*numHazards-1)+numPHParams+1+numHazards)+1:(2^Mval-1)+(2^Mval-1)*(hazCtr-1))], 
					paste(outfilename, '/', namesHazGroups[hazCtr], sep = ''), betanames = betaNames)
			row.names(outputdata$summary)[1:(2^Mval+1)] = 
					paste(row.names(outputdata$summary)[1:(2^Mval+1)], namesHazGroups[1], sep = '_')
			row.names(outputdata$d) = paste(row.names(outputdata$d), namesHazGroups[1], sep = '_')
			row.names(outputdata$H) = paste('H', namesHazGroups[1], sep = '_')
			row.names(outputdata$Rmp) = paste(row.names(outputdata$Rmp), namesHazGroups[1], sep = '_')
		} else {
			tempoutputdata = AnalyzeComplete(Mval, censortime, 
				finaldata[,c(1:2^Mval+2^Mval*(hazCtr-1) + 1,
					2^Mval*(2*numHazards-1)+numPHParams+1+hazCtr, 
					(2^Mval*(2*numHazards-1)+numPHParams+1+numHazards)+1:(2^Mval-1)+(2^Mval-1)*(hazCtr-1))], 
					paste(outfilename, '/', namesHazGroups[hazCtr], sep = ''))
			row.names(tempoutputdata$summary)[1:(2^Mval+1)] = 
				paste(row.names(tempoutputdata$summary)[1:(2^Mval+1)], namesHazGroups[hazCtr], sep = '_')
			row.names(tempoutputdata$d) = paste(row.names(tempoutputdata$d), namesHazGroups[hazCtr], sep = '_')
			row.names(tempoutputdata$H) = paste('H', namesHazGroups[hazCtr], sep = '_')
			row.names(tempoutputdata$Rmp) = paste(row.names(tempoutputdata$Rmp), namesHazGroups[hazCtr], sep = '_')
			outputdata$summary = rbind(outputdata$summary[1:(hazCtr-1),], tempoutputdata$summary[1,],
									   outputdata$summary[1:(2^Mval*(hazCtr-1))+hazCtr-1,], 
									   tempoutputdata$summary[1:2^Mval+1,],
									   outputdata$summary[(2^Mval*(hazCtr-1)+hazCtr):nrow(outputdata$summary),])
			outputdata$d = rbind(outputdata$d, tempoutputdata$d)
			outputdata$H = rbind(outputdata$H, tempoutputdata$H)
			outputdata$Rmp = rbind(outputdata$Rmp, tempoutputdata$Rmp)
		}
	}
	# Calculation of AIC, BIC, DIC and add likelihood to chains
	nrow.data = nrow(finaldata)
	binwidth = censortime/2^Mval
	Rmpvec = as.data.frame(matrix(NA, ncol = (2^Mval-1)*2*numHazards, nrow = nrow.data))
	Rmpvec[,1:((2^Mval-1)*2*numHazards/2)*2-1] = 
			finaldata[,1:((2^Mval-1)*numHazards) + 1+2^Mval*numHazards+numPHParams+2^Mval*(numHazards-1)+numHazards]
	Rmpvec[,1:((2^Mval-1)*2*numHazards/2)*2] = 
			1-finaldata[,1:((2^Mval-1)*numHazards) + 1+2^Mval*numHazards+numPHParams+2^Mval*(numHazards-1)+numHazards]
	logds = NULL
	for(hazctr in 1:numHazards){
		logds = cbind(logds, 
					  matrix(matrix(rep(log(finaldata[,1+numPHParams+2^Mval*(numHazards*2-1)+hazctr]), 
										each = 2^Mval), ncol = 2^Mval, nrow = nrow.data, byrow = TRUE)+ 
							 t(mat01%*%t(log(Rmpvec[,1:((2^Mval-1)*2)+(2^Mval-1)*2*(hazctr-1)]))) - 
							 log(binwidth), nrow = nrow.data))
	}
	if(numPHParams > 0){	
		Xbeta = rowSums(as.matrix(Xfixed[rep(1:n, nrow.data),])*(as.matrix(finaldata[,1:numPHParams+2^Mval*numHazards+1])[rep(1:nrow.data, each = n),]))
	} else { Xbeta = rep(0, n*nrow.data)	}
	hazBin = failBin[,rep(1:2^Mval, numHazards)]
	cumulBin = inBin[,rep(1:2^Mval, numHazards)]
	for(hazctr in 1:numHazards){
		hazBin[(1:n)[-indices[[hazctr]]], 1:2^Mval + (hazctr-1)*2^Mval] = 0
		cumulBin[(1:n)[-indices[[hazctr]]], 1:2^Mval + (hazctr-1)*2^Mval] = 0
	}
	totalone = rep(delta, nrow.data)*(rowSums(hazBin[rep(1:n, nrow.data),]*(logds[rep(1:nrow.data, each = n),]))+Xbeta) -
					rowSums(cumulBin[rep(1:n, nrow.data),]*(finaldata[,1:(2^Mval*numHazards)+1][rep(1:nrow.data, each = n),]))*exp(Xbeta)
	logliks = rowSums(matrix(totalone, ncol = n, byrow = TRUE))
	finaldata = as.data.frame(cbind(finaldata, logliks))
	names(finaldata) = c(names(finaldata)[-ncol(finaldata)], 'loglik')
	write.table(finaldata, paste(outfilename, '/MCMCchains.txt', sep = ''), row.names = FALSE)
	# The number of estimated parameters is (alpha, lambda)*number of haz + number of Rmps estimated + H00*number of haz + number PH covs
	numberofEstParams = 2*numHazards + numRmpsEstimated + 2*numHazards + numPHParams
	AIC = 2*(numberofEstParams) - 2*min(logliks)
	BIC = -2*min(logliks) + (numberofEstParams)*log(n)
	DIC = .5*var(-2*logliks) + mean(-2*logliks)
	
	hazardRate = outputdata$d/(censortime/2^Mval)
	colnames(hazardRate) = paste('q', c('.5', '.025', '.975'), sep = '')
	rownames(hazardRate) = paste(rep(paste('h.Bin', 1:2^Mval, sep = ''), numHazards), 
								 rep(paste('Group', 1:numHazards, sep = ''), each = 2^Mval), sep = '.')

	if(GR == FALSE){ gr.used = 'Option not used' } else { gr.used = 'Option Used'	}
	outputdata = c(outputdata, list(hazardRate = hazardRate, AIC = AIC, BIC = BIC, DIC = DIC, 
									burnIn = burnIn, thin = thin, TotalIters = loopctr-1, 
									convergence = convergence, gelman.rubin.used = gr.used, fix.thin = fix.thin,
									fix.burnIn = fix.burnIn, fix.max = fix.max, initialValues = initialValues,
									runtime = paste(round((unclass(Sys.time()) - unclass(systemtime))/(60*60), 2), 'hours'), 
									outfolder = outfilename, maxStudyTime = censortime))
	if(numPHParams > 0){
		GraphNPbetas(M = Mval, dests = outputdata$d, np.betaests = outputdata$beta[-(1:numPHParams),], 
					 numhazgroups = numHazards, hazgroupnames = namesHazGroups, censortime = censortime, 
					 file = paste(outfilename, '/', sep = ''))
	} else {
		GraphNPbetas(M = Mval, dests = outputdata$d, np.betaests = outputdata$beta, 
					 numhazgroups = numHazards, hazgroupnames = namesHazGroups, censortime = censortime, 
					 file = paste(outfilename, '/', sep = ''))
	}
	class(outputdata) = "MRH"
	return(outputdata)
			
}

