#' @title Apply sine model on the full set of genes or isoforms 
#' @usage OscopeSine(DataInSc, parallel=FALSE, parallelParam=NULL)
#' @param DataInSc a gene-by-sample (isoform-by-sample) matrix indicating the rescaled expression of two genes/isoforms.
#' all values should be bettwen [-1, 1].
#' @param parallel whether apply parallel computing. if it is TRUE, BiocParallel will be called.
#' @param parallelParam a SnowParam object to specify the clusters. If it is NULL, the default
#' will be set as SnowParam(workers = 5, type = "SOCK")
#' remains unchanged for over NCThre iterations.
#' @return Output is a list with 4 sublists, each shows a N-by-N matrix, in which
#' N is the total number of genes (isoforms).
#' SimiMat: similarity matrix (sine scores); the sine scores are calculated
#' by -log10(epsilon^2).
#' DiffMat: dissimilarity matrix; shown are epsilon^2 for each gene pair.
#' ShiftMat: optimal phase shift estimate for each pair of genes.
#' @examples aa <- sin(seq(0,1,.1))
#' bb <- sin(seq(0.5,1.5,.1))
#' cc <- sin(seq(0.9,1.9,.1))
#' OscopeSine(rbind(aa,bb,cc))
#' @author Ning Leng


OscopeSine <- function(DataInSc, parallel=FALSE,parallelParam=NULL){
	expect_is(DataInSc, "matrix")
	if(is.null(rownames(DataInSc)))stop("No gene/isoform names!")
	if(length(unique(rownames(DataInSc)))!=nrow(DataInSc)) stop("Duplicated gene/isoform names!")

	NumGene <- nrow(DataInSc)

	if(parallel==FALSE)
		Res <- sapply(1:(NumGene-1),function(i)SineFun(DataInSc, i),simplify=FALSE)

  if(parallel){
		if(is.null(parallelParam))parallelParam <- SnowParam(workers = 5, type = "SOCK")
		Res <- bplapply(1:(NumGene-1),function(i)Oscope::SineFun(DataInSc, i),BPPARAM=parallelParam)
	}
		expect_is(Res, "list")
		Out <- FormatSineOut(Res, DataInSc)[1:3]
}



