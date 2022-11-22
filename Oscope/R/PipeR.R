#' @title Calculate residual of polynomial fit
#' @usage PipeR(Data,Ndg=3,Method="Poly")
#' @param Data gene-by-sample matrix or isoform-by-sample matrix.It should be rescaled to values bwteen
#' [-1,1].
#' @param Ndg degree of polynomial.
#' @param Method only polynomial fitting ("Poly") is available now.
#' @return The function will fit polynomial curve to each row of the data. 
#' The output returns the MSE of each row (gene/isoform).
#' @examples aa <- sin(seq(0,1,.1))
#' bb <- sin(seq(0.5,1.5,.1))
#' cc <- sin(seq(0.9,1.9,.1))
#' res <- PipeR(rbind(aa,bb,cc))
#' @author Ning Leng

PipeR <- function(Data,Ndg=3,Method="Poly"){
	V <- 1:ncol(Data)
	# MSE
	if(Method=="Poly")ResV <- t(sapply(1:nrow(Data),function(i)mean(residuals(lm(Data[i,]~poly(V,Ndg)))^2)))
	colnames(ResV) <- rownames(Data)
	ResV
}

