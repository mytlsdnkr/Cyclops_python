#' @title Function for searching optimal phase shift
#' @usage SineOptim(Pairdata)
#' @param Pairdata a sample-by-2 matrix indicating the rescaled expression of two genes/isoforms.
#' all values should be bettwen [-1, 1].
#' @return Output provides the optimal phi (shift) and its corresponding epsilon^2 (value) of the sine model. 
#' epsilon_{g1,g2}^2 = sum_s [X_{g1,s}^2+X^2_{g2,s}
#'                           -2X_{g1,s}X_{g2,s} cos(phi_{g1,g2})-sin^2(phi_{g1,g2})]^2
#' @examples aa <- sin(seq(0,1,.1))
#' bb <- sin(seq(0.5,1.5,.1))
#' SineOptim(cbind(aa,bb))
#' @author Ning Leng
SineOptim <- function(Pairdata){
	out <- try(optim(par=0, 
									 function(Pairdata,par){sum(abs(Pairdata[,2]^2 + Pairdata[,1]^2 - 2 * Pairdata[,1] * Pairdata[,2] * cos(par) - (sin(par))^2))}
									 ,Pairdata=Pairdata, lower=0,upper=2*pi,method="L-BFGS-B"))
	V <- c(NA,NA)
	if(class(out)!="try-error"){
		V <- c(out$value,out$par)
	}
	#expect_is(V, c("integer","numeric"))
	names(V) <- c("value","shift")
	V
}

