#' @title Run k-medoid algorithm with varying k on similarity matrix
#' @usage scanK(SimiMatIn, quan=.95,cut=NULL, maxK=NULL,minSize=0, maxSize=200, fixK=NULL, rawscale=FALSE)
#' @param SimiMatIn gene-by-gene similarity matrix
#' @param quan only gene pairs with similarity score >= quan th quantile will be 
#' considered in the cluster analyses. Default is 0.95.
#' @param cut pre-defined cutoff. Gene pairs with similarity score >= cut will be considered in
#' cluster analyses. If cut is defined, quan will be ignored.
#' @param maxK max number of clusters to consider (scan). if numbC=NULL, it will be calculated as 
#' [number of gene considered]/10.
#' @param minSize,maxSize Only clusters with minSize<= cluster size <= maxSize are 
#' reported in output.
#' @param fixK if fixK is specified, the k-medoids algorithm will be applied with fixK clusters.
#' @param rawscale 
#' Recall the input 
#' is the similarity matrix (-log10(distance from the sine model)). 
#' the k-medoids clustering will be applied using (-Input) as distance. If rawscale is defined as TRUE,
#' the k-medoids clustering will be applied using -10^Input as distance.
#' @return scanK() function runs k-medoid clustering with varying number of clusters (k). 
#' The k is varied from 2 to maxK. The input of scanK() function should be a similarity matrix.
#' scanK() function will cluster genes in gene pairs with high similarity score (the threshold can be
#' defined using parameter quan). To select the top genes, the function first calculate the max similarity
#' score for each gene, then select the genes with high max score.
#'
#' The output object is a list with 4 sublists:
#' membOut: members in each cluster. clusters are sorted by median similarity score within cluster;
#'
#' MedCor: median similarity score for each cluster;
#'
#' Mat: input similarity matrix;
#'
#' filteredMat: similarity matrix, only showing the top genes used in clustering;
#'
#' Kcluster: cluster indicator of each top gene.
#' @examples aa <- sin(seq(0,1,.1))
#' bb <- sin(seq(0.5,1.5,.1))
#' cc <- sin(seq(0.9,1.9,.1))
#' tmp <- matrix(sin(rnorm(330)),ncol=11)
#' rownames(tmp) <- paste0("tmp",1:30)
#' Dat <- rbind(aa, bb, cc, tmp)
#' res1 <- OscopeSine(Dat)
#' res2 <- scanK(res1$SimiMat, quan=.8, maxK=5)
#' @author Ning Leng

scanK <- function(SimiMatIn, quan=.95, cut=NULL, maxK=NULL,minSize=0, maxSize=200,fixK=NULL,
rawscale=FALSE){
	if(is.null(rownames(SimiMatIn))) stop("Row names are not provided!")
	if(is.null(colnames(SimiMatIn))) stop("Column names are not provided!")
	if(length(unique(rownames(SimiMatIn)))!=nrow(SimiMatIn)) stop("Duplicated gene/isoform names!")	
	expect_is(SimiMatIn, "matrix")
	expect_equal(nrow(SimiMatIn),ncol(SimiMatIn))
	#library(cluster)
	
	# RM ones with Inf
	WhichInf <- which(rowSums(abs(SimiMatIn))==Inf)
	if(length(WhichInf)>0)SimiMatIn <- SimiMatIn[-WhichInf,-WhichInf]
	cor_max <- sapply(1:nrow(SimiMatIn),function(i)max(SimiMatIn[i,-i],na.rm=TRUE))
	QQ <- quantile(cor_max,quan,na.rm=TRUE)
	if(!is.null(cut))QQ <- cut
	expect_is(QQ,c("numeric","integer"))
	message("gene pairs above this threshold are considered:")
	message(QQ)
	wcm <- which(cor_max>=QQ)
	sMatCor <- SimiMatIn[wcm, wcm]
	expect_is(sMatCor, "matrix")

	numC <- maxK
	if(is.null(fixK)){
	if (is.null(numC))
		    numC <- ceiling((dim(SimiMatIn)[1]/10)*(1-quan))
	if(numC<2) numC <- 2
	  message("max number of clusters considered:",numC)
	
	

	TryList <- vector("list",numC-1)
	asw <- rep(NA,numC-1)
	for(i in 2:numC){
	if(rawscale==FALSE){
	TryList[[i-1]] <- pam(dist(-sMatCor), k=i)}
	if(rawscale==TRUE){
	temp <- 10^(-sMatCor)
	diag(temp)=0
	TryList[[i-1]] <- pam(dist(temp), k=i)}
	asw[i-1] <- TryList[[i-1]]$ silinfo $ avg.width
	}
	k.best <- which.max(asw)
	expect_is(k.best, "integer")

	message("optimal number of clusters:", k.best+1)
	Try <- TryList[[k.best]]
	}
	if(!is.null(fixK)){ 
	if(rawscale==FALSE)Try <- pam(dist(-sMatCor), k=fixK )
	if(rawscale==TRUE){
	temp <- 10^(-sMatCor)
        diag(temp) <- 0
	Try <- pam(dist(temp), k=fixK )}
}
	memb <- Try$clustering

	a <- table(memb)
	a2 <- a[which(a>=minSize & a <=maxSize)]
	membOut <- sapply(1:length(a2),function(i)names(memb)[which(memb==names(a2)[i])],simplify=FALSE)

	MeanCor <- sapply(membOut,function(i)median(SimiMatIn[i,i]))
	expect_is(MeanCor,c("numeric","integer"))
	Order <- order(MeanCor,decreasing=TRUE)
	membOutSort <- membOut[Order]
	MeanCorSort <- MeanCor[Order]
	outNames <- paste0("cluster",1:length(MeanCor))
	names(membOutSort) <- outNames
	names(MeanCorSort) <- outNames
	return(Out <- list(membOut=membOutSort,MedCor=MeanCorSort,
		Mat=SimiMatIn,filteredMat=sMatCor,
		Kcluster=memb))
}
