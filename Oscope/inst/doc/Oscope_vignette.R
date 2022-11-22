### R code from vignette source 'Oscope_vignette.Rnw'

###################################################
### code chunk number 1: style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: Oscope_vignette.Rnw:61-62
###################################################
library(Oscope)


###################################################
### code chunk number 3: Oscope_vignette.Rnw:88-91
###################################################
data(OscopeExampleData)
str(OscopeExampleData)
set.seed(10)


###################################################
### code chunk number 4: Oscope_vignette.Rnw:101-102
###################################################
Sizes <- MedianNorm(OscopeExampleData)


###################################################
### code chunk number 5: Oscope_vignette.Rnw:110-111
###################################################
DataNorm <- GetNormalizedMat(OscopeExampleData, Sizes)


###################################################
### code chunk number 6: Oscope_vignette.Rnw:130-133
###################################################
MV <- CalcMV(Data = OscopeExampleData, Sizes = Sizes)
str(MV$GeneToUse)
DataSubset <- DataNorm[MV$GeneToUse,]


###################################################
### code chunk number 7: Oscope_vignette.Rnw:153-156 (eval = FALSE)
###################################################
## MV2 <- CalcMV(Data = DataNorm, Sizes = NULL, NormData = TRUE)
## str(MV2$GeneToUse)
## DataSubset2 <- DataNorm[MV2$GeneToUse,]


###################################################
### code chunk number 8: Oscope_vignette.Rnw:174-175
###################################################
DataInput <- NormForSine(DataNorm)


###################################################
### code chunk number 9: Oscope_vignette.Rnw:195-197 (eval = FALSE)
###################################################
## SineRes <- OscopeSine(DataInput)
## str(SineRes)


###################################################
### code chunk number 10: Oscope_vignette.Rnw:203-205
###################################################
SineRes <- OscopeSine(DataInput, parallel=TRUE)
str(SineRes)


###################################################
### code chunk number 11: Oscope_vignette.Rnw:217-219 (eval = FALSE)
###################################################
## DataInput2 <- NormForSine(DataSubset)
## SineRes2 <- OscopeSine(DataInput2)


###################################################
### code chunk number 12: Oscope_vignette.Rnw:227-229
###################################################
KMRes <- OscopeKM(SineRes, maxK = 10)
print(KMRes)


###################################################
### code chunk number 13: Oscope_vignette.Rnw:268-269
###################################################
ToRM <- FlagCluster(SineRes, KMRes, DataInput)


###################################################
### code chunk number 14: Oscope_vignette.Rnw:271-274
###################################################
print(ToRM$FlagID_bysine)
print(ToRM$FlagID_byphase)
print(ToRM$FlagID) # all flagged clusters 


###################################################
### code chunk number 15: Oscope_vignette.Rnw:283-285
###################################################
KMResUse <- KMRes[-ToRM$FlagID]
print(KMResUse)


###################################################
### code chunk number 16: Oscope_vignette.Rnw:294-296 (eval = FALSE)
###################################################
## ENIRes <- OscopeENI(KMRes = KMResUse, Data = DataInput, NCThre = 100)
## print(ENIRes)


###################################################
### code chunk number 17: Oscope_vignette.Rnw:302-303 (eval = FALSE)
###################################################
## ENIRes <- OscopeENI(KMRes = KMResUse, Data = DataInput, NCThre = 100, parallel=TRUE)


###################################################
### code chunk number 18: Oscope_vignette.Rnw:320-321 (eval = FALSE)
###################################################
## DataNorm2 <- DataNorm[,ENIRes[["cluster2"]]]


###################################################
### code chunk number 19: Oscope_vignette.Rnw:331-336 (eval = FALSE)
###################################################
## par(mfrow = c(3,2))
## for(i in 1:6)
## plot(DataNorm[KMResUse[["cluster2"]][i], ENIRes[["cluster2"]]],
## xlab = "Recovered order", ylab = "Expression",
## main = KMResUse[["cluster2"]][i])


###################################################
### code chunk number 20: Oscope_vignette.Rnw:348-353 (eval = FALSE)
###################################################
## par(mfrow = c(3,2))
## for(i in 1:6)
## plot(DataNorm[KMResUse[["cluster3"]][i], ENIRes[["cluster3"]]],
## xlab = "Recovered order", ylab = "Expression",
## main = KMResUse[["cluster3"]][i])


###################################################
### code chunk number 21: Oscope_vignette.Rnw:361-362
###################################################
print(sessionInfo())


