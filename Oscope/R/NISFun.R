#' @title Run Extended Nearest Insertion and 2-opt on a gene cluster identified by OscopeKM function
#' @usage NISFun(ClusterList, DataIn, i, Ndg=3, NChun=4, RdmStart=FALSE, N=20000, NCThre=1000)
#' @param DataIn gene-by-sample matrix or isoform-by-sample matrix.It should be rescaled to values bwteen
#' [-1,1].
#' @param i the cluster of interest. If the second cluster in ClusterList is of interest, specify i=2.
#' @param ClusterList a list of gene clusters. Each sublist contains a vector of gene names.
#' @param Ndg degree of polynomial.
#' @param NChun number of starting points for polynomial fitting.
#' @param RdmStart whether the start points are randomly selected.
#' @param N,NCThre The 2-opt algorithm will stop if N iterations has been performed or if the optimal order 
#' remains unchanged for over NCThre iterations.
#' @return This function performs the extended nearest insertion (ENI) and 2-opt algorithm to
#' a particular cluster identified by OscopeKM function.
#' The ENI algorithm searchs for the optimal sample order which minimizes the MSE of 
#' sliding polynomial regression (SPR). 
#' This function will call PipeShiftCDF() function, which fits 
#' SPR to each row of the data. 
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
#' The 2-opt algorithm is then applied to improve the optimal order searching of the ENI.
#' In each iteration, 2-opt algorithm will randomly choose two points (samples), the flip the samples
#' between these two points. The new order will be adapted if it provides smaller SPR MSE.
#' The output returns the optimal order for the cluster of interest. 
#' @examples aa <- sin(seq(0,1,.1))
#' bb <- sin(seq(0.5,1.5,.1))
#' cc <- sin(seq(0.9,1.9,.1))
#' res <- NISFun(list(c("aa","bb"),"cc"), rbind(aa,bb,cc),i=1, NChun=2, N=50)
#' @author Ning Leng


NISFun <- function(ClusterList, DataIn,i, Ndg=3, NChun=4, RdmStart=FALSE, N=20000, NCThre=1000){
	 expect_is(ClusterList,"list")
	 UnkIn <- ClusterList[[i]]
	 NamesTmp <- UnkIn
	 expect_that(str(NamesTmp), prints_text("chr"))
	 if(length(setdiff(NamesTmp, rownames(DataIn)))>0)stop("Some genes in the cluster list don't have a corresponding expression entry!")
	 message("ENI of ",names(ClusterList)[i], "; inserting samples:")
	 Outtmp <- ImpShift(DataIn[NamesTmp,],Ndg=Ndg, NChun=NChun, RdmStart=RdmStart)
	 message("2-opt of ",names(ClusterList)[i])
	 OutList <-Opt2Shift(DataIn[NamesTmp,], N=N, NCThre=NCThre, Seq=Outtmp)
	 Out <- OutList[[1]]	   
	 expect_is(Out, c("integer","numeric"))
	 Out
			   }

