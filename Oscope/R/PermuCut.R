#' @title Define sine scroe cutoff using permuted data
#' @usage PermuCut(Data, NumPermu=1000)
#' @param Data a gene-by-sample (isoform-by-sample) matrix indicating the rescaled expression of genes/isoforms.
#' all values should be between [-1, 1].
#' @param NumPermu number of permuted genes to generate.
#' @return Output contains a vector of numbers. Each number presents max sine score of a given 
#' permuted gene.
#' @examples aa <- sin(seq(0,1,.1))
#' bb <- sin(seq(0.5,1.5,.1))
#' cc <- sin(seq(0.9,1.9,.1))
#' tmp <- matrix(sin(rnorm(330)),ncol=11)
#' rownames(tmp) <- paste0("tmp",1:30)
#' Dat <- rbind(aa, bb, cc, tmp)
#' res1 <- PermuCut(Dat,100)
#' @author Ning Leng


PermuCut <- function(Data, NumPermu=1000){
	RP <- TRUE
	if(NumPermu<=nrow(Data))RP <- FALSE
	S <- sample(1:nrow(Data),NumPermu,replace=RP)
	Ncol <- ncol(Data)
	Data1 <- t(sapply(1:NumPermu, function(i)sample(Data[S[i],],Ncol)))
	rownames(Data1) <- paste0("permu",1:NumPermu)

	SineRes1 <- OscopeSine(Data1)
	SimiMatIn <- SineRes1$SimiMat

	cor_max <- sapply(1:nrow(SimiMatIn),function(i)max(SimiMatIn[i,-i],na.rm=TRUE))
	out <- list(SimiMat=SimiMatIn,MaxEach=cor_max)	
}
