#' @title	Search for the optimal sample order by using the Extended Nearest Insertion 
#' @usage ImpShift(Data, Seq=NULL, NChun=4, RdmStart=FALSE, Ndg=3)
#' @param Data gene-by-sample matrix or isoform-by-sample matrix.It should be rescaled to values bwteen
#' [-1,1].
#' @param Ndg degree of polynomial.
#' @param NChun number of starting points for polynomial fitting.
#' @param RdmStart whether the start points are randomly selected.
#' @param Seq NULL or a vector indicates the sample order.
#' if specified, the samples will be first reordered by this vector.
#' @return This function performs the extended nearest insertion (ENI).
#' The ENI algorithm searchs for the optimal sample order which minimizes the MSE of 
#' sliding polynomial regression  (SPR). 
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
#' The output returns the optimal order which provides the smallest SPR MSE.
#' @examples aa <- sin(seq(0,1,.1))
#' bb <- sin(seq(0.5,1.5,.1))
#' cc <- sin(seq(0.9,1.9,.1))
#' res <- ImpShift(rbind(aa,bb,cc), NChun=2)
#' @author Ning Leng
ImpShift <- function(Data, Seq=NULL, NChun=4, RdmStart=FALSE, Ndg=3){
	expect_is(Data, "matrix")
	if(is.null(Seq))Seq=1:ncol(Data)	
	expect_is(Seq, "integer")
	Od <- Seq[1:3]
	Ncol <- ncol(Data)
	pb <- txtProgressBar(min=4, max=Ncol, style=3)
	for(i in 4:Ncol){
	tp <- Seq[i]
	#matrix with i columns, i rows
 		
	mat <- t(sapply(1:(i-2),function(j)c(Od[1:j],tp,Od[(j+1):(i-1)])))
	mat <- rbind(c(tp, Od[1:(i-1)]),mat,c( Od[1:(i-1)], tp))
	#if(is.null(Ndg))Ndg <- ceiling(i/Seg)
	# 11/16 error when Ndg > 27
	if(Ndg>27)Ndg <- 27
	# or input Ndg
	Stat <- sapply(1:i,function(j)PipeShiftCDF(Data[,mat[j,]], Ndg=Ndg, RdmStart=RdmStart, NChun=NChun))
	Min <- which.min(Stat)
	Od <- mat[Min,]
	setTxtProgressBar(pb, i)
	}
	close(pb)
Od}

