#' @title Rescale the gene/isoform expression matrix 
#' @usage NormForSine(Data, qt1=.05, qt2=.95)
#' @param Data input gene-by-sample matrix or isoform-by-sample matrix
#' @param qt1,qt2 thresholds for outlier adjustment. For each gene/isoform,
#' values <= qt1 th quantile (>= qt2 th quantile)
#' will be pushed to qt1 th quantile (qt2 th quantile) prior to the scaling.
#' default values are 0.05 and 0.95.
#' @return The output will be a gene-by-sample or isoform-by-sample matrix.
#' For each gene/isoform, the expressions will be scaled linearly to [-1,1]
#' @examples NormForSine(matrix(rnorm(10), nrow=2))
#' @author Ning Leng


NormForSine <- function(Data, qt1=.05, qt2=.95){
	expect_is(Data,"matrix")	
	Q5 <- apply(Data,1,function(i)quantile(i,.05))
	Q95 <- apply(Data,1,function(i)quantile(i,.95))
	Rg <- Q95-Q5
	expect_is(Rg,c("numeric","integer"))
	DataSc2 <- ((Data-Q5)*2/Rg)-1
	DataSc2[which(DataSc2<(-1),arr.ind=TRUE)] <- -1
	DataSc2[which(DataSc2>1,arr.ind=TRUE)] <- 1
	expect_is(DataSc2,"matrix")
	Out <- DataSc2
}
