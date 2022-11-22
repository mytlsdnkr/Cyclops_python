#' @title Calculate residual of the sliding polynomial regression
#' @usage PipeShiftCDF(Data, Ndg=3, NChun=4, RdmStart=FALSE)
#' @param Data gene-by-sample matrix or isoform-by-sample matrix. It should be rescaled to values bwteen
#' [-1,1].
#' @param Ndg degree of polynomial.
#' @param NChun number of starting points for polynomial fitting.
#' @param RdmStart whether the start points are randomly selected.
#' @return The function will fit sliding polynomial regression  (SPR)
#' to each row of the data. 
#' For each gene/isoform, SPR fits NChun polynomial curves with different starting 
#' points (samples). The samples with smaller order than the start point will be appended
#' to follow the last sample when fitting. So each fitting consider same number
#'  of samples. If RdmStart = TRUE, the start points are randomly selected.
#' Otherwise they are evenly sampled along the sample order.
#' The aggregated MSE of a fit (using a specific start point) is defined as the 
#' summation of the MSEs of all genes/isoforms considered here.
#' The output returns the MSE of the SPR, which is the largest aggregated MSE across fits 
#' using different start points.
#' @examples aa <- sin(seq(0,1,.1))
#' bb <- sin(seq(0.5,1.5,.1))
#' cc <- sin(seq(0.9,1.9,.1))
#' res <- PipeShiftCDF(rbind(aa,bb,cc), NChun=2)
#' @author Ning Leng
PipeShiftCDF <- function(Data,Ndg=3, NChun=4, RdmStart=FALSE){
	#expect_is(Data, "matrix")
	N <- ncol(Data)
	Start <- seq(from=1,to=N,by=floor(N/NChun))
	if(RdmStart==TRUE & N>=4){
		Num <- sample(0:(Start[2]-2),1)
		Start[-c(1,length(Start))] <- Start[-c(1,length(Start))]-Num
	}
	Seq0 <- sapply(2:NChun,function(i)c(Start[i]:N, 1:(Start[i]-1)))
	Seq <- rbind(1:N, t(Seq0))
	t0 <- sapply(1:NChun,function(i)PipeR(Data[,Seq[i,]], Ndg),simplify= FALSE)
  EC <- sapply(1:NChun, function(i)quantile(ecdf(t0[[i]]), 1:100/100))
	out0 <- colSums(EC)
	out <- max(out0)
	expect_is(out,c("integer","numeric"))
	out
}

