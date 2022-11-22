#' @title Calculate absolute correlations among gene pairs
#' @usage AbsCor(DataIn, method="pearson", diagNA=TRUE)
#' @param DataIn input data, gene-by-sample matrix
#' @param method "pearson" or "spearman"; default is "pearson"
#' @param diagNA whether replace diagonal values to NA's
#' @return Output is a gene-by-gene matrix; 
#' the i, j th entry shows the absolute correlation of the ith and jth gene.
#' @examples AbsCor(matrix(rnorm(10),ncol=5))
#' @author Ning Leng
AbsCor <- function(DataIn, method = "pearson", diagNA=TRUE){
		MatCor <- cor(scale(t(DataIn)))
    absMatCor <- abs(MatCor)
		diag(absMatCor) <- NA
		SimiMatIn <- absMatCor
}
