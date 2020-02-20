#!/usr/bin/env Rscript

##-------------------------------------------------------------------
# Author: Nilgun Donmez
# Date: December 22, 2015
# Last modified by Salem Malikic in May 2019
#
# Part of CTPsingle package. See README.md for more information.
##-------------------------------------------------------------------

library(DPpackage)
library(lpSolve)


usage <- function()
{
	warning("Usage: CTPsingle_pancancer.R -f <mutation_file> -o <output_prefix> -m <adjacency_matrices> [ -r <num_restarts> -a <alpha> ]")
	quit()
}

getAncestry <- function(parents, numNodes)
{
	ancestry <- matrix(0, numNodes, numNodes)	
	for(i in 1:numNodes)
	{
		ancestry[i,i] <- 1		# every node is an ancestor of itself (by definition)
		j <- parents[i]
		while(j>0)				# zero means it has no ancestor (i.e. root node)
		{
			ancestry[j,i] <- 1
			j <- parents[j]
		}
	}
	return(ancestry)
}

assignClusters <- function(clones, clusters)
{
	nClones <- length(clones)
	nClusters <- length(clusters)
	costs <- matrix(0, nClusters, nClones)
	
	for(i in 1:nClusters)	
		costs[i,] <- abs(clones - clusters[i])

	fit <- lp.assign(costs)	
	return(fit)
}

assignBetas <- function(solutions, clusters, numNodes, children)
{
	nVars <- 2*numNodes
	coeffs <- rep(0, nVars)		# this will be the objective function
	coeffs[(numNodes+1):nVars] <- 1		# these are the delta variables
	
	nConstraints <- 1 + 3*numNodes
	constraints <- matrix(0, nConstraints, nVars)		
	dirs <- rep(">=", nConstraints)
	rhs <- rep(0, nConstraints)
	
	constraints[1, 1] <- 1		# this is the root node
	rhs[1] <- 1.0
	dirs[1] <- "<="				# the root node frequency shall not exceed 1.0
	
	j <- 2
	for(i in 1:numNodes)
	{
		constraints[j, 1:numNodes] <- solutions[i,]				# delta_i + beta_i >= Freq_i
		constraints[j, (numNodes+1):nVars] <- solutions[i,]
		
		constraints[j+1, 1:numNodes] <- (0 - solutions[i,])		# delta_i - beta_i >= -Freq_i
		constraints[j+1, (numNodes+1):nVars] <- solutions[i,]
				
		rhs[j] <- clusters[i]
		rhs[j+1] <- (0.0 - clusters[i])
		
		constraints[j+2, 1:numNodes] <- (0 - children[i,])		# beta_parent - beta_child1 - beta_child2 >= 0
		constraints[j+2, i] <- 1		
		rhs[j+2] <- 0
		
		j <- j + 3
	}		
	fit <- lp(direction="min", objective.in=coeffs, const.mat=constraints, const.dir=dirs, const.rhs=rhs)
}

applyBootstrap <- function(clusters, ancestry, children, numNodes, restarts, tolerance=0.001)	
{	
	results <- list()
	results$score <- 2 * numNodes		# the objective can not be more than this
	for(k in 1:restarts)
	{
		alphas <- rep(1.0/numNodes, numNodes)	# always try this first
		if(k>1)	
		{
			alphas <- runif(numNodes, 0.0, 0.1)
			alphas <- alphas/sum(alphas)
		}				
		betas <- rep(0, numNodes)
		for(i in 1:numNodes)	
			betas[i] <- sum(ancestry[i,] * alphas)	
		
		mappings <- matrix(0, numNodes, numNodes)
		prevObj <- 2 * numNodes		
		converged <- FALSE
		while(converged == FALSE)
		{
			fit <- assignClusters(betas, clusters)			
			cat("Step1: ", fit$objval, "\n")		
			mappings <- fit$solution
			
			fit2 <- assignBetas(mappings, clusters, numNodes, children)
			cat("Step2: ", fit2$objval,"\n")
			betas <- fit2$solution[1:numNodes]		# the rest is deltas' which we don't need
			
			if(tolerance > abs(prevObj-fit2$objval))
			{	
				cat("Converged...\n")
				converged <- TRUE			
				if(fit2$objval < results$score)
				{
					results$mappings <- mappings
					results$betas <- betas
					results$score <- fit2$objval					
				}
			}
			prevObj <- fit2$objval					
		}
	}
	return(results)
}

solveTreeILP <- function(clusters, adjMatFile, outputPrefix, numNodes, restarts)
{
	x <- as.matrix(read.table(adjMatFile, header=F, skip=2))
	x <- x + 1		# so that the indices start from 1
	
	for(i in 1:nrow(x))
	{
		parents <- rep(0, numNodes)
		children <- matrix(0, numNodes, numNodes)
		j <- 2
		while(j < ncol(x))
		{
			parents[x[i,j+1]] <- x[i,j]
			children[x[i,j], x[i,j+1]] <- 1
			j <- j + 2
		}		
		ancestry <- getAncestry(parents, numNodes)	
		ilpSolution <- applyBootstrap(clusters, ancestry, children, numNodes, restarts)	
	
		bestSolution <- cbind(ilpSolution$mappings, signif(clusters, digits=2), parents, 1:numNodes, signif(ilpSolution$betas, digits=2), rep(ilpSolution$score, numNodes))  	
		write.table(bestSolution, paste(outputPrefix, ".treeSize_", numNodes, "-treeIndex_", i, ".CTPsingle.tree", sep=""), quote=F, row.names=F, col.names=F, sep=" ")
	}	
}

dpCluster <- function(inputFile, outputPrefix, adjFolder, alpha, a0, b0, a1, b1, ngrid, nburn, nsave, restarts, useAlpha=FALSE, maxMutations=10010, minMutations=50, maxDownsampling=0.11, maxFreq=1.05, minMembers_parameter=20)
{
	frqData <- read.table(inputFile, header=TRUE)	
	numMuts = nrow(frqData)
	minMembers = 0.03*numMuts
	if(minMembers < minMembers_parameter)
	{
		minMembers = minMembers_parameter
	}

	if(numMuts < minMutations)
	{	
		warning("No or too few mutations found. Exiting without clustering.")
		cat(file=paste(outputPrefix, ".CTPsingle.purity", sep=""), "NA")
		quit()
	}	
	if(numMuts > maxMutations)
	{
		if((maxMutations/numMuts) < maxDownsampling)
		{
			warning("Hypermutated dataset. Exiting without clustering.")
			quit()
		}
		warning("Downsampling the dataset.")
		frqData <- frqData[sample(1:numMuts, maxMutations, replace=FALSE),]		
		numMuts <- nrow(frqData)
	}	
	
	frqData[,9] <- rep(2, numMuts)
	if(frqData[1,8] == "male")
	{
		frqData[frqData[1,] == "X", 9] <- 1
		frqData[frqData[1,] == "Y", 9] <- 1
	}
	
	y <- as.matrix(frqData[,5:6])										# mutant count, reference count	
	y[,2] <- ((2/frqData[,7]) * y[,1]) + ((2/frqData[,9]) * y[,2])		# total count after adjusting based on copy number of the mutation
	

	if(numMuts > nburn)		
	{
		nburn <- numMuts + 10
		nsave <- floor(1.2*nburn)
		warning("Increasing nburn and nsave to match the size of input data.")
	}

	priorList = list(a0=a0, b0=b0, a1=a1, b1=b1)
	if(useAlpha)
		priorList = list(alpha=alpha, a1=a1, b1=b1)
	dpbb <- DPbetabinom(y, ngrid=ngrid, prior=priorList, mcmc=list(nburn=nburn, nsave=nsave, nskip=3, ndisplay=100), state=NULL, status=T)

	numClusters <- dpbb$state$ncluster
	Assignments <- dpbb$state$ss

	num <- 0
	mostLikely <- Assignments
	for(i in 1:numClusters)	
	{
		mu <- mean((2*y[Assignments==i,1])/y[Assignments==i,2])	
		total <- sum(Assignments==i)		
		if(mu < maxFreq && total > minMembers)		
		{
			num <- num + 1
			mostLikely[Assignments==i] <- num		
		}
		else
			mostLikely[Assignments==i] <- 0		# this is a ghost cluster
	}
	numClusters <- num		# adjust to the new number
	Assignments <- mostLikely	

	clusters <- rep(0, numClusters)	
	
	for(i in 1:numClusters)	
		clusters[i] <- mean((2*y[Assignments==i,1])/y[Assignments==i,2])			

	tumorPurity <- max(clusters)
	cat(file=paste(outputPrefix, ".CTPsingle.purity", sep=""), "Purity: \t", tumorPurity)
					
	clusters <- clusters/tumorPurity		# adjust cluster frequencies by tumor purity	
	averageFrequency <- rep(0.0, numMuts)
	for(i in 1:numClusters)	
		averageFrequency[mostLikely==i] <- clusters[i]
			
	totalCN <- frqData[,7]
	mutationCN <- rep(1, numMuts)
	certainty <- cbind(frqData[,1:4], mostLikely, averageFrequency, totalCN, mutationCN)				
	write.table(certainty, paste(outputPrefix, ".CTPsingle.clusters", sep=""), quote=F, row.names=F, sep="\t")

	if(numClusters > 1 && numClusters < 10)		
	{
		adjMatFile <- paste(adjFolder,"/AdjacencyMatrix", numClusters, ".txt", sep="")
		solveTreeILP(clusters, adjMatFile, outputPrefix, numClusters, restarts)
	}
	return(tumorPurity)
}


processCmds <- function(myargs)
{	
	inputFile <- ""
	outputPrefix <- ""
	
	# default parameters
	alpha <- 0.001
	a0 <- 1.2
	b0 <- 0.01
	a1 <- 5.0
	b1 <- 5.0
	ngrid <- 100
	nburn <- 1000
	nsave <- 2000
	restarts <- 100
	useAlpha <- TRUE
	adjFolder <- "../GammaAdjMatrices"

	cmnd <- 1
	while(cmnd < length(myargs))
	{
		if(myargs[cmnd][[1]] == "-f")
			inputFile <- myargs[cmnd+1][[1]]
			
		if(myargs[cmnd][[1]] == "-o")
			outputPrefix <- myargs[cmnd+1][[1]]
			
		if(myargs[cmnd][[1]] == "-a")
			alpha <- as.numeric(myargs[cmnd+1][[1]])	
			
		if(myargs[cmnd][[1]] == "-a0")
		{
			a0 <- as.numeric(myargs[cmnd+1][[1]])
			useAlpha <- FALSE
		}
			
		if(myargs[cmnd][[1]] == "-b0")
		{
			b0 <- as.numeric(myargs[cmnd+1][[1]])
			useAlpha <- FALSE
		}
			
		if(myargs[cmnd][[1]] == "-a1")
			a1 <- as.numeric(myargs[cmnd+1][[1]])
			
		if(myargs[cmnd][[1]] == "-b1")
			b1 <- as.numeric(myargs[cmnd+1][[1]])
			
		if(myargs[cmnd][[1]] == "-r")
			restarts <- as.numeric(myargs[cmnd+1][[1]])
			
		if(myargs[cmnd][[1]] == "-g")
			ngrid <- as.numeric(myargs[cmnd+1][[1]])
			
		if(myargs[cmnd][[1]] == "-b")
			nburn <- as.numeric(myargs[cmnd+1][[1]])

		if(myargs[cmnd][[1]] == "-s")
			nsave <- as.numeric(myargs[cmnd+1][[1]])
			
		if(myargs[cmnd][[1]] == "-m")
			adjFolder <- myargs[cmnd+1][[1]]
			
		cmnd <- cmnd + 2
	}
	
	tumorPurity <- dpCluster(inputFile, outputPrefix, adjFolder, alpha, a0, b0, a1, b1, ngrid, nburn, nsave, restarts, useAlpha)
	cat("Estimated Purity: ", tumorPurity, "\n")	
}

argLine <- commandArgs(TRUE)
myargs <- (strsplit(argLine, " "))

if(length(myargs) < 3)
	usage()

processCmds(myargs)
proc.time()
warnings()


