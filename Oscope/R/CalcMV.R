#' @title Calculate estimated mean and variance of RNA-Seq data
#' @usage CalcMV(Data, Sizes=NULL, NormData=FALSE, MeanCutLow=100, MeanCutHigh=NULL, ApproxVal=10^-6, Plot=TRUE)
#' @param Data input data matrix; it should be a gene-by-sample or isoform-by sample matrix
#' @param Sizes The library size factor for each sample. the number of values in Sizes is expected to be the
#' same as the number of columns of Data. The library size factor will be estimated using the median
#' normalization method implemented in EBSeq if Sizes is specified as NULL.
#' @param MeanCutLow,MeanCutHigh we suggests the users to apply Oscope on genes with high mean and 
#' high variance. By default, MeanCutLow is specified as 100, consequently only genes with mean > 100
#' will be used. The CalcMV function will fit a linear regression on log(variance)~log(mean) on these 
#' genes. Genes with variance above this line are considered as the high mean high variance genes.
#' The upper bound of mean may be specified using MeanCutHigh. If both are specified as NULL, all of the genes
#' will be considered when fitting the regression.
#' @param NormData whether the data is already normalized. If NormData=TRUE, the specification of Sizes
#' will be ignored and no normalization will be applied.
#' @param ApproxVal Default is 10^-6. It is used to approximate the estimate of parameter q for genes/isoforms
#' whose estimated variance is less than estimated mean. q will be estimated using 1-ApproxVal
#' @param Plot if Plot = T, a mean-variance plot will be shown. The fitted line will be shown and the
#' selected genes will be marked in green.
#' @return Output is a list with 6 sublists : Mean: estimated means of genes/isoforms; Var: estimated variances;
#' Median: estimated medians; GeneToUse: the high mean high variance genes (suggested input for Oscope);
#' Q: estimated q's (without apporximation); Q_mdf: estimated q's with approximations;
#' Phi_mdf: estimated overdispersion parameter (phi), with approximations.
#' @examples 
#' exp=matrix(rnorm(100,1000,10),ncol=10)
#' rownames(exp)=paste0("g",1:10)
#' CalcMV(exp)
#' @author Ning Leng

CalcMV<-function(Data, Sizes=NULL, NormData=FALSE, MeanCutLow=100, MeanCutHigh=NULL, ApproxVal=10^-6, Plot=TRUE){
	expect_is(Data, "matrix")
	expect_is(rownames(Data), "character")
	EmpData <- Data
	EmpSizes <- Sizes
	if(is.null(Sizes))EmpSizes <- MedianNorm(EmpData)
	expect_is(EmpSizes,c("numeric","integer"))

	# Calculate normalized means
	if(NormData==FALSE)EmpData.norm <- t(t( EmpData )/EmpSizes)  
	else EmpData.norm <- EmpData
	MeansC1 <- rowMeans(EmpData.norm)
	MedC1 <- apply(EmpData.norm,1,median)
	expect_is(MeansC1,c("numeric","integer"))
	expect_is(MedC1,c("numeric","integer"))

	# Calculate var
	Sig_tmp <- (EmpData-MeansC1%*%t(EmpSizes))^2
	Sig_tmp2 <- t(t(Sig_tmp)/EmpSizes)
	VarC1 <- rowMeans(Sig_tmp2)
	expect_is(VarC1,c("numeric","integer"))

	# calculate q 
	QC1 <- MeansC1/VarC1
	
	# Some genes are with mean>=var (q>=1)
	# In this case, use 1-10^-6 to approximate q
	QNB <- QC1
	QNB[which(QNB>=1)] <- ApproxVal
	
	# calculate phi
	PhiNB <- (1-QNB)/(MeansC1*QNB)
	PhiInput <- PhiNB
	# option to simulate constant phi for all genes

	MVOut <- list(Mean=MeansC1,Var=VarC1, Median=MedC1, 
						Q=QC1, Q_mdf=QNB, Phi_mdf=PhiNB)
	SampleMean <- MVOut$Mean
	SampleVar <- MVOut$Var

	Which <- 1:length(SampleMean)
	if(!is.null(MeanCutLow))
		{
			Which <- which(SampleMean>MeanCutLow)
			if(!is.null(MeanCutHigh))
				Which <- intersect(Which,which(SampleMean<MeanCutHigh))
		}
	if(length(Which)<3)stop("Too few genes are selected based on the settings of MeanCutHigh and MeanCutLow!")	
	Meanfit <- log10(SampleMean[Which])
	Varfit <- log10(SampleVar[Which])
	lm1 <- lm(Varfit~Meanfit)
	expect_is(lm1, "lm")
	Coef <- coef(lm1)

	Gt10 <- names(SampleMean)[Which]
	MeanUse <- SampleMean[Gt10]
	VarUse <- SampleVar[Gt10]
	Fit <- 10^(Coef[1] + Coef[2]*log10(MeanUse))
	expect_is(Fit, c("numeric","integer"))
	Diff <- VarUse-Fit
	SamplePickGenes <- names(MeanUse)[which(Diff>0)]

	if(Plot==TRUE){
	plot(SampleMean, SampleVar, col="gray",pch=21, xlab="Mean",
			 ylab="Variance", log="xy")
	lines(MeanUse[order(MeanUse)], Fit[order(MeanUse)])
	points(MeanUse[SamplePickGenes], VarUse[SamplePickGenes], pch=21,col="green")
	abline(v=c(MeanCutLow, MeanCutHigh))}

	out<-list(Mean=MeansC1,Var=VarC1, Median=MedC1, GeneToUse=SamplePickGenes,  
						            Q=QC1, Q_mdf=QNB, Phi_mdf=PhiNB)

}
