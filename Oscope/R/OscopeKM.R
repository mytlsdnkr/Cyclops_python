#' @title Oscope K medoid module
#' @usage OscopeKM(SineRes, quan=.95,cut=NULL,maxK=NULL,minSize=0, maxSize=200, fixK=NULL, rawscale=TRUE)
#' @param SineRes output of OscopeSine function. 
#' @param quan only gene pairs with similarity score >= quan th quantile will be 
#' considered in the clustering analyses. Default is 0.95.
#' @param cut pre-defined cutoff. Gene pairs with similarity score >= cut will be considered in
#' cluster analyses. If cut is defined, quan will be ignored.
#' @param maxK max number of clusters to consider (scan). if numbC=NULL, it will be calculated as 
#' [number of gene considered]/10
#' @param minSize,maxSize Only clusters with minSize<= cluster size <= maxSize are 
#' reported in output.
#' @param fixK if fixK is specified, the k-medoids algorithm will be applied with fixK clusters.
#' @param rawscale 
#' Recall the input 
#' is the similarity matrix (-log10(distance from the sine model)). 
#' the k-medoids clustering will be applied using (-Input) as distance. If rawscale is defined as TRUE,
#' the k-medoids clustering will be applied using -10^Input as distance.
#' @return OscopeKM() calls scanK() function, which runs k-medoid clustering with varying number of clusters (k). 
#' The k is varied from 2 to maxK. The input should be the output of OscopeSine() function.
#' scanK() function will cluster genes in gene pairs with high similarity score (the threshold can be
#' defined using parameter quan). To select the top genes, the function first calculate the max similarity
#' score for each gene, then select the genes with high max score.
#'
#' The output object shows
#' members in each cluster. clusters are sorted by median similarity score within cluster.
#' @examples aa <- sin(seq(0,1,.1))
#' bb <- sin(seq(0.5,1.5,.1))
#' cc <- sin(seq(0.9,1.9,.1))
#' tmp <- matrix(sin(rnorm(330)),ncol=11)
#' rownames(tmp) <- paste0("tmp",1:30)
#' Dat <- rbind(aa, bb, cc, tmp)
#' res1 <- OscopeSine(Dat)
#' res2 <- OscopeKM(res1, quan=.8, maxK=5)
#' @author Ning Leng


OscopeKM <- function(SineRes, quan=.95, cut=NULL, maxK=NULL,minSize=0, maxSize=200,fixK=NULL,rawscale=TRUE){
	InMat <- SineRes$SimiMat
	expect_is(InMat, "matrix")
	expect_equal(nrow(InMat),ncol(InMat))
	membSine <- scanK(SimiMatIn=InMat, quan=quan, cut=cut, maxK=maxK,
										minSize=minSize, maxSize=maxSize,fixK=fixK,rawscale=rawscale)
	Out <- membSine[[1]]
}
