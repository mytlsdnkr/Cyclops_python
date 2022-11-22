#' @title Run the 2-opt algorithm to improve the optimal order searching of the Extended Nearest Insertion
#' @usage Opt2Shift(Data,N=20000,Seq,Ndg=3,NChun=4, NCThre=1000, RdmStart=FALSE)
#' @param Data gene-by-sample matrix or isoform-by-sample matrix. It should be rescaled to values bwteen
#' [-1,1].
#' @param N,NCThre The 2-opt algorithm will stop if N iterations has been performed or if the optimal order 
#' remains unchanged for over NCThre iterations.
#' @param Ndg degree of polynomial.
#' @param NChun number of starting points for polynomial fitting.
#' @param RdmStart whether the start points are randomly selected.
#' @param Seq a vector indicates the sample order obtained from the ENI.
#' @return This function performs the 
#' the 2-opt algorithm to improve the optimal order searching of the Extended Nearest Insertion (ENI).
#' In each iteration, the function will randomly choose two points (samples), the flip the samples
#' between these two points. The new order will be adapted if it provides smaller SPR MSE.
#' The output returns the optimal order and its SPR MSE. 
#' @examples aa <- sin(seq(0,1,.1))
#' bb <- sin(seq(0.5,1.5,.1))
#' cc <- sin(seq(0.9,1.9,.1))
#' res <- ImpShift(rbind(aa,bb,cc), NChun=2)
#' res2 <- Opt2Shift(rbind(aa,bb,cc), NChun=2, N=50, Seq=res)
#' @author Ning Leng




Opt2Shift <- function(Data,N=20000,Seq,Ndg=3,NChun=4, NCThre=1000, RdmStart=FALSE){
	expect_is(Data, "matrix")
	Ncol <- ncol(Data)
	expect_is(Seq, "integer")
	SeqIter <- Seq
	StatIter <- PipeShiftCDF(Data[,Seq],Ndg=Ndg,RdmStart=RdmStart)
	nc <- 0
	for(i in 1:N){
		Choose <- sample(1:Ncol,2)
		Seq0 <- SeqIter
		Seq0[Choose[1]:Choose[2]] <- SeqIter[Choose[2]:Choose[1]]
		Stat <- PipeShiftCDF(Data[,Seq0],Ndg=Ndg,RdmStart=RdmStart, NChun=NChun)
		nc <- nc+1
		if(Stat<StatIter){
			  SeqIter <- Seq0
		  	StatIter <- Stat
			  nc <- 0
			  message("updated at iteration " ,i, " ; ", appendLF=FALSE)
		}
		if(nc>NCThre){
			message("last iteration:",i)
			break
		}}
	expect_is(SeqIter,"integer")
	expect_is(StatIter,c("integer","numeric"))
Out <- list(SeqIter, StatIter)
}



