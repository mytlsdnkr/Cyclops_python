#' @title Format SinFun outputs from lists to matrix
#' @usage FormatSineOut(result, DataInSc, ShiftRg=pi/4)
#' @param result Output from SineFun
#' @param DataInSc a gene-by-sample (isoform-by-sample) matrix indicating the rescaled expression of two genes/isoforms.
#' all values should be bettwen [-1, 1].
#' @param ShiftRg phase shift cutoff. 
#' @return Output is a list with 4 sublists, each shows a N-by-N matrix, in which#' N is the total number of genes (isoforms).
#' SimiMat: similarity matrix (sine scores); the sine scores are calculated
#' by -log10(epsilon^2).
#' DiffMat: dissimilarity matrix; shown are epsilon^2 for each gene pair.
#' ShiftMat: optimal phase shift estimate for each pair of genes.
#' @examples aa <- sin(seq(0,1,.1))
#' bb <- sin(seq(0.5,1.5,.1))
#' cc <- sin(seq(0.9,1.9,.1))
#' DataInSc <- rbind(aa,bb,cc)
#' NumGene <- nrow(DataInSc)
#' Res <- sapply(1:(NumGene-1),function(i)SineFun(DataInSc, i),simplify=FALSE)
#' Out <- FormatSineOut(Res, DataInSc)
#' @author Ning Leng

FormatSineOut <- function(result, DataInSc, ShiftRg=pi/4){
	expect_is(DataInSc, "matrix")
	expect_is(result, "list")
	DataIn <- DataInSc
	if(length(result)<(length(result[[1]][[1]])+1))
		result <- c(result,list(NULL))# sometimes one short while using sapply
	NumGene <- length(result)

	MatKS <- matrix(0,ncol=NumGene,nrow=NumGene)
	for(i in 1:(NumGene-1)){
	MatKS[i,(i+1):NumGene] <- -log10(result[[i]]$value)
	MatKS[(i+1):NumGene,i] <- -log10(result[[i]]$value)	
	}
	MatDiff <- 10^(-MatKS)
	diag(MatDiff) <- 0
	MatKSShift <- matrix(0,ncol=NumGene,nrow=NumGene)
	for(i in 1:(NumGene-1)){
	MatKSShift[i,(i+1):NumGene] <- result[[i]]$shift
	MatKSShift[(i+1):NumGene,i] <- result[[i]]$shift
	}

	rownames(MatKS)=colnames(MatKS)=rownames(MatDiff)=colnames(MatDiff)=
	rownames(MatKSShift)=colnames(MatKSShift)=rownames(DataIn)
	# push the non-shift pairs with score NA
 	LeftPi <- MatKSShift%%pi
	WhichNA <- which(LeftPi < ShiftRg | LeftPi > pi-ShiftRg, arr.ind=TRUE)
	Mat_cut <- MatKS
	if(length(WhichNA)>0)Mat_cut[WhichNA] <- NA	
	Out <- list(SimiMat=MatKS, DiffMat=MatDiff,
		 ShiftMat=MatKSShift,SimiMat_ShiftOnly=Mat_cut)
}
