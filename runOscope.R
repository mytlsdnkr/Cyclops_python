####use Oscope to do selection on the eigen genes and output the eigen clusters used for CYCLOPS
###select groups of possible periodic eigenenes using Oscope
library(Oscope)
library(BiocParallel)
args <- commandArgs(trailingOnly = T)

print(args)

indir <- "./Result/"
dir(indir)
infiles <- grep("EigengeneExp.csv", dir(indir, full.names = TRUE), value = TRUE)
print(infiles)
outfiles <- gsub("EigengeneExp.csv", "EigengeneExp_OscopeCluster.csv", infiles)
print(outfiles)

FormatSineOut <- function(result, DataInSc, ShiftRg=pi/4){
# 	expect_is(DataInSc, "matrix")
# 	expect_is(result, "list")
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
SineOptim <- function(Pairdata){
	out <- try(optim(par=0,	function(Pairdata,par){sum(abs(Pairdata[,2]^2 + Pairdata[,1]^2 - 2 * Pairdata[,1] * Pairdata[,2] * cos(par) - (sin(par))^2))}
									 ,Pairdata=Pairdata, lower=0,upper=2*pi,method="L-BFGS-B"))
	V <- c(NA,NA)
	if(class(out)!="try-error"){
		V <- c(out$value,out$par)
	}
	#expect_is(V, c("integer","numeric"))
	names(V) <- c("value","shift")
	V
}

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


OscopeSine <- function(DataInSc, parallel=FALSE,parallelParam=NULL){
	if(is.null(rownames(DataInSc)))stop("No gene/isoform names!")
	if(length(unique(rownames(DataInSc)))!=nrow(DataInSc)) stop("Duplicated gene/isoform names!")

	NumGene <- nrow(DataInSc)

	if(parallel==FALSE)
		Res <- sapply(1:(NumGene-1),function(i)SineFun(DataInSc, i),simplify=FALSE)

  if(parallel){
		if(is.null(parallelParam))parallelParam <- SnowParam(workers = 5, type = "SOCK")
		Res <- bplapply(1:(NumGene-1),function(i)Oscope::SineFun(DataInSc, i),BPPARAM=parallelParam)
	}
# 		expect_is(Res, "list")
		Out <- FormatSineOut(Res, DataInSc)[1:3]
}



NormForSine <- function(Data, qt1=.05, qt2=.95){
	Q5 <- apply(Data,1,function(i)quantile(i,.05))
	Q95 <- apply(Data,1,function(i)quantile(i,.95))
	Rg <- Q95-Q5
	DataSc2 <- ((Data-Q5)*2/Rg)-1
	DataSc2[which(DataSc2<(-1),arr.ind=TRUE)] <- -1
	DataSc2[which(DataSc2>1,arr.ind=TRUE)] <- 1
	Out <- DataSc2
}
ii <- 1
for (ii in 1:length(infiles))  {
  infile <- infiles[ii]
  outfile <- outfiles[ii]
  eigenD <- read.csv(infile)
  OscopeData <- as.matrix(eigenD[,-c(1,2)])

  dimnames(OscopeData) <- list("r"=eigenD[,2], "c"=colnames(eigenD)[-c(1,2)] )
  DataInput <- NormForSine(OscopeData)
  ##apply sine model on the eigengenes
  SineRes <- OscopeSine(DataInput, parallel = TRUE)

  SineRes
#   KMRes <- try(OscopeKM(SineRes, quan=0.4, maxK = 10), silent = TRUE)     
#   outD <- NULL
#   class(KMRes)
#   if (class(KMRes) == "list")  {
#     ##check the number of each element
#     cnum = unlist(lapply(KMRes, length))
#     cindex = as.numeric(which(cnum > 1))
#     if (length(cindex) < length(cnum))  {
#       KMRes = KMRes[cindex]
#     }
#     ##Flag clusters with small within-cluster sine scores and/or small within-cluster phase shifts
#     ToRM <- FlagCluster(SineRes, KMRes, DataInput)
#     KMResUse <- KMRes[-ToRM$FlagID]
#     if (length(KMResUse))  {
#       eigenName <- sapply(KMResUse, function(z) {paste(z, collapse = "|")} )
#       eigenIndex <- sapply(KMResUse, function(z) { gz <- gsub("eigen_(\\d+)_\\S+", "\\1", z, perl = TRUE)
#       return(paste(sort(gz), collapse = "|"))
#       } )
#       outD <- data.frame(groupName = names(KMResUse), eigenName = eigenName, eigenIndex = eigenIndex)
#     }
#   }
#   ##select those pair eigen genes with defined cutoff
#   ##it may need to adjust this cutoff according to specific data
#   cutoff <- as.numeric(quantile(SineRes$DiffMat, probs = seq(0, 1, by=0.1) )[2] )
#   print(quantile(SineRes$DiffMat, probs = seq(0, 1, by=0.1) ))
#   similarD <- SineRes$DiffMat
#   peigenName <- peigenIndex <- NULL
#   eigenID <- rownames(similarD)
#   rown <- length(eigenID)
#   for (i in 1:(rown-1) )  {
#     for (j in (i+1):(rown) )  {
#       if (similarD[j,i] <= cutoff)  {
#         peigenName <- c(peigenName , paste(eigenID[i], eigenID[j], sep = "|") )
#         peigenIndex <- c(peigenIndex, paste(sort(c(i,j)), collapse ="|"))
#       }
#     }
#   }
#   ##get the selected pair cluster
#   pairD <- NULL
#   if (length(peigenName))  {
#     pairD <- data.frame(groupName = paste("cluster0", 1:length(peigenName), sep=""),
#                         eigenName = peigenName, eigenIndex = peigenIndex )
#   }
#   ##output the Oscope selected and pair cluster
#   outD <- rbind(outD, pairD)
#   if (length(outD))  {
#     ##get the duplicated rows
#     dup <- which(duplicated(outD$eigenIndex) == TRUE)
#     if (length(dup))  {
#       write.csv(outD[-dup,], file = outfile, row.names = FALSE)
#     }  else  {
#       write.csv(outD, file = outfile, row.names = FALSE)
#     }
#   }
}
