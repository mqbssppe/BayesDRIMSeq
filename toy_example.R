source('BayesDRIMSeq.R')

#	We consider a toy example with 100 genes

nGenes <- 100
genes <- c()
for(i in 1:nGenes){
	nIso <- 1 + rpois(n=1,  lambda = 2)		# number of transcripts per gene
	genes <- c(genes, rep(paste0("GENE_",i), nIso))
}
genes <- as.factor(genes)

nTranscripts <- length(genes)				# total number of transcripts
nSamplesA <- 3						# number of replicates for 1st condition 
nSamplesB <- 3						# number of replicates for 2nd condition
nSamples <- nSamplesA + nSamplesB			# total number of replicates
# 	Simulate data.frame with counts (no DTU evidence at all)
countDataFrame <- data.frame( matrix(rpois(n = nTranscripts*nSamples, lambda = 30), nrow = nTranscripts, ncol = nSamples) )

# 	Call the BayesDRIMSEQ function as:
myRes <- laplaceDM(
		count_data = countDataFrame, 		# data.frame of counts with dimension: nTranscripts x nSamples
		gene_data = genes, 			# factor with `nGenes` levels with length: nTranscripts
		grouping = as.factor(
			c(rep('A',nSamplesA),
				rep('B',nSamplesB))), 	# factor with 2 levels and length nSamples
		min_reads_filter = 20, 			# positive integer used to filter our low expressed transcripts
		nCores = 8, 				# number of paraller workers
		lambdaRate = 0.5			# positive prior parameter \lambda
	)


# 	RESULTS: Our FDR decision rule d_4 corresponds to `myRes$fdrTrust`
