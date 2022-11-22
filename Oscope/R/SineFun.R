#' @title Apply sine model on one particular gene vs. other genes
#' @usage SineFun(DataInSc,i)
#' @param DataInSc a gene-by-sample (isoform-by-sample) matrix indicating the rescaled expression of two genes/isoforms.
#' all values should be bettwen [-1, 1].
#' @param i the gene (isoform) of interest. The function will apply the sine model on gene (isoform) i vs. 
#' gene (isoform) j for all j > i. Gene (isoform) i (j) is defined as the gene (isoform )shown in the i (j) th
#' row. i should be smaller than the total number of genes (isoforms).
#' @return Output is a list with two sublists, each shows the optimal phi's (shift) and epsilon's (value).
#' N-i entries will be included in each sublist (N is the total number of genes/isoforms). The kth entry
#' indicates results of gene (isoform) i vs. i+k.
#' @examples aa <- sin(seq(0,1,.1))
#' bb <- sin(seq(0.5,1.5,.1))
#' cc <- sin(seq(0.9,1.9,.1))
#' SineFun(rbind(aa,bb,cc), 1)
#' @author Ning Leng

SineFun <- function(DataInSc,i){
		#expect_is(DataInSc, "matrix")
		#expect_is(i, c("numeric","integer"))
		out <- sapply((i+1):nrow(DataInSc),function(j){
		a1 <- SineOptim(cbind(DataInSc[i,],DataInSc[j,]))
		a1
			},simplify=FALSE)
		  out0 <- list(value=sapply(out,function(k)k[1]),
		            shift=sapply(out,function(k)k[2]))
		  }

