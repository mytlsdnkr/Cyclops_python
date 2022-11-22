#' @title Search for the optimal sample order for different gene clusters
#' @usage OscopeENI(KMRes, Data, ClusterUse=NULL, Ndg=3, NChun=4, RdmStart=FALSE,
#' N=20000, NCThre=1000, parallel=FALSE, parallelParam=NULL)
#' @param Data gene-by-sample matrix or isoform-by-sample matrix. It should be rescaled to values bwteen
#' [-1,1].
#' @param KMRes output of OscopeKM() function.
#' @param ClusterUse a vector indicating what clusters are of interest. For example, by setting ClusterUse
#' = c(1,2,3), only the top 3 clusters will be considered while recovering the base cycle order. 
#' If ClusterUse=NULL, all clusters will be used.
#' @param Ndg degree of polynomial.
#' @param NChun number of starting points for polynomial fitting.
#' @param RdmStart whether the start points are randomly selected.
#' @param N,NCThre The 2-opt algorithm will stop if N iterations has been performed or if the optimal order 
#' @param parallel whether apply parallel computing. if it is TRUE, BiocParallel will be called.
#' @param parallelParam a SnowParam object to specify the clusters. If it is NULL, the default
#' will be set as SnowParam(workers = 5, type = "SOCK")
#' remains unchanged for over NCThre iterations.
#' @return This function performs the extended nearest insertion (ENI) and 2-opt algorithm to
#' all clusters (or a subset of picked clusters) identified by OscopeKM function.
#' The function will recover independent orders to each of the clusters.
#' For each cluster, the ENI algorithm will be applied to
#' search for the optimal sample order which minimizes the MSE of 
#' sliding polynomial regression (SPR). 
#' This function will call PipeShiftCDF() function, which fits 
#' SPR to expression of each gene/isoform within a cluster. 
#' For each gene/isoform, SPR fits NChun polynomial curves with different starting 
#' points (samples). The samples with smaller order than the start point will be appended
#' to follow the last sample when fitting. So each fitting consider same number
#'  of samples. If RdmStart = TRUE, the start points are randomly selected.
#' Otherwise they are evenly sampled along the sample order.
#' The aggregated MSE of a fit (using a specific start point) is defined as the 
#' summation of the MSEs of all genes/isoforms considered here.
#' The MSE of the SPR is defined as the largest aggregated MSE across fits 
#' using different start points.
#' The output of PipeShiftCDF() returns the optimal order which provides the smallest SPR MSE.
#' The 2-opt algorithm will then be applied to improve the optimal order searching of the ENI.
#' In each iteration, the 2-opt algorithm will randomly choose two points (samples), the flip the samples
#' between these two points. The new order will be adapted if it provides smaller SPR MSE.
#' The output returns the optimal order for each cluster of interest. 
#' It is a list with multiple sublists, in which each sublist includes the recovered order
#' of the corresponding cluster in ClusterUse. If ClusterUse is not specified, the k th sublist
#' shows the recovered order in KMRes
#' @examples aa <- sin(seq(0,1,.1))
#' bb <- sin(seq(0.5,1.5,.1))
#' cc <- sin(seq(0.9,1.9,.1))
#' dd <- sin(seq(1.2,2.2,.1))
#' res <- OscopeENI(list(c1=c("aa","bb"),c2=c("cc","dd")), rbind(aa,bb,cc,dd), NChun=2, N=50)
#' @author Ning Leng



OscopeENI <- function (KMRes, Data, ClusterUse=NULL, Ndg=3, NChun=4, RdmStart=FALSE, N=20000, NCThre=1000, parallel=FALSE, parallelParam=NULL){
	if(!is.list(KMRes))stop("Input KMRes is not a list!")
	if(is.null(names(KMRes)))stop("Please specify cluster names (list names) ")
	if(is.null(ClusterUse)) ClusterUse <- 1:length(KMRes)
	expect_is(Data, "matrix")
	if(parallel==FALSE)
	Res <- sapply(1:length(ClusterUse), function(k)NISFun(KMRes, Data, i=ClusterUse[k], Ndg=Ndg,
																					NChun=NChun, RdmStart=RdmStart, N=N, NCThre=NCThre ),
						 															simplify=FALSE)
	
	if(parallel){
	if(is.null(parallelParam))parallelParam <- SnowParam(workers = 5, type = "SOCK")
	expect_is(parallelParam, "SnowParam")
	Res <- bplapply(1:length(ClusterUse), function(k)Oscope::NISFun(KMRes, Data, i=ClusterUse[k], Ndg=Ndg,
																					NChun=NChun, RdmStart=RdmStart, N=N, NCThre=NCThre ),
						 															BPPARAM=parallelParam)
	
	}
	expect_is(Res, "list")
	expect_equal(length(KMRes),length(Res))
	names(Res) <- names(KMRes)[ClusterUse]
	Res
}
