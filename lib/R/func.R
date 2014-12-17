#!/usr/bin/Rscript
#args <- commandArgs(T)
#qcFile <- args[1]
#output <- args[2]

#read.delim(qcFile, header=T) -> qc

addQual <- function(qc) {

qc$del <- as.character(as.vector(qc$del))
qc$ins <- as.vector(qc$ins)
qc$pVal <- 1
qc$pValAdj <- 1
#qc$pValPois <- 1
qc$qScore <- 0
qc$qScoreI <- 0
qc$qScoreD <- 0

pVal.I <- rep(1, nrow(qc))
pVal.D <- rep(1, nrow(qc))

delSel <- qc$del != "."
data.frame(matrix(unlist(strsplit(unlist(strsplit(as.vector(qc$del[delSel]), "\\|")), ":")), ncol=2, byrow=T)) -> dels
dels$ln <- nchar(as.vector(dels$X1))
colnames(dels) <- c("seq", "occ", "ln")
delPees <- ppois(1:50, mean(rep(dels$ln, dels$occ)), lower.tail=F)

insSel <- qc$ins != "."
data.frame(matrix(unlist(strsplit(unlist(strsplit(as.vector(qc$ins[insSel]), "\\|")), ":")), ncol=2, byrow=T)) -> inss
inss$ln <- nchar(as.vector(inss$X1))
colnames(inss) <- c("seq", "occ", "ln")
insPees <- ppois(1:50, mean(rep(dels$ln, dels$occ)), lower.tail=F)

# Calculate errorrate per amplicon
amps <- unique(qc$amp)
for (i in amps) {
	#print(i)
	sel <- qc$amp == i
        prob <- sum(qc$nVar[sel]) / (sum(qc$nRef[sel]) + sum(qc$nVar[sel]))
	#prob <- median(qc$nVar[sel] / (qc$nRef[sel] + qc$nVar[sel]))
        pees <- pbinom(qc$nVar[sel], qc$nRef[sel] + qc$nVar[sel], prob, lower.tail=F)
        qc$pVal[sel] <- pees
#       pois <- ppois(qc$nVar[sel], (qc$nRef + qc$nVar) * prob, lower.tail=F)
#       qc$pValPois[sel] <- pois

	# Seperate for indel/deletion events
	selDels <-  delSel & sel
	if (sum(selDels) > 0) {
		data.frame(matrix(unlist(strsplit(unlist(strsplit(as.vector(qc$del[selDels]), "\\|")), ":")), ncol=2, byrow=T)) -> dels2
		dels2$ln <- nchar(as.vector(dels2$X1))
		colnames(dels2) <- c("seq", "occ", "ln")
		dels2$occ <- as.numeric(as.vector(dels2$occ))
	}
	for (j in which(selDels)) {
		data.frame(matrix(unlist(strsplit(unlist(strsplit(as.vector(qc$del[j]), "\\|")), ":")), ncol=2, byrow=T)) -> dels
		dels$ln <- nchar(as.vector(dels$X1))
		colnames(dels) <- c("seq", "occ", "ln")
		dels$occ <- as.numeric(as.vector(dels$occ))
		#dels$pDel <- delPees[as.integer(dels$ln)]
		#dels$pVal <- ppois(as.numeric(dels$occ), dels$pDel, lower.tail=F)
		dels$pDel <- sapply(dels$ln, function(x) {sum(dels2$occ[dels2$ln == x]) / (sum(qc$nRef[sel]) + sum(qc$nVar[sel])) })
		dels$pVal <- pbinom(dels$occ, qc$dp[j], dels$pDel, lower.tail=F)
		dels$Q <- -10 * log10(dels$pVal)
			
		high <- order(dels$occ, -dels$Q, decreasing=T)[1]
		qc$del[j] <- paste(dels$seq[high], dels$occ[high], sep=":")		
		qc$qScoreD[j] <- dels$Q[high]
		pVal.D[j] <- dels$pVal[high]
	}
	selIns <-  insSel & sel
	if (sum(selIns) > 0) {
		data.frame(matrix(unlist(strsplit(unlist(strsplit(as.vector(qc$ins[selIns]), "\\|")), ":")), ncol=2, byrow=T)) -> inss2
		inss2$ln <- nchar(as.vector(inss2$X1))
		colnames(inss2) <- c("seq", "occ", "ln")
		inss2$occ <- as.numeric(as.vector(inss2$occ))
	}
	for (j in which(selIns)) {
		data.frame(matrix(unlist(strsplit(unlist(strsplit(as.vector(qc$ins[j]), "\\|")), ":")), ncol=2, byrow=T)) -> inss
		inss$ln <- nchar(as.vector(inss$X1))
		colnames(inss) <- c("seq", "occ", "ln")
		inss$occ <- as.numeric(as.vector(inss$occ))
		#inss$pIns <- insPees[as.integer(inss$ln)]
		#inss$pVal <- ppois(as.numeric(inss$occ), inss$pIns, lower.tail=F)
		inss$pIns <- sapply(inss$ln, function(x) {sum(inss2$occ[inss2$ln == x]) / (sum(qc$nRef[sel]) + sum(qc$nVar[sel])) })
		inss$pVal <- pbinom(inss$occ, qc$dp[j], inss$pIns, lower.tail=F)
		inss$Q <- -10 * log10(inss$pVal)

		high <- order(inss$occ, -inss$Q, decreasing=T)[1]
		qc$ins[j] <- paste(inss$seq[high], inss$occ[high], sep=":")
		qc$qScoreI[j] <- inss$Q[high]
		pVal.I[j] <- inss$pVal[high]
	}
	#print(sum(-10 * log10(pees) > 100))
	#print(sum(qc$qScoreI[selIns] > 100))
	#print(sum(qc$qScoreD[selDels] > 100))
}

p.adjust(qc$pVal) -> qc$pValAdj
-10 * log(p.adjust(pVal.I)) -> qc$qScoreI
-10 * log(p.adjust(pVal.D)) -> qc$qScoreD
qc$qScore <- -10 * log10(qc$pValAdj)
#qc$qScore[qc$pValAdj == 0] <- 1000
qc$qScore[qc$qScore > 1000] <- 1000
qc$qScoreI[qc$qScoreI > 1000] <- 1000
qc$qScoreD[qc$qScoreD > 1000] <- 1000
qsel <- (qc$qScore >= 100 & qc$nVar >= 10) | (qc$qScoreI >= 200) | (qc$qScoreD >= 200)
dpfilt <- as.numeric(qc$dp) >= 100
#varfilt <- as.numeric(qc$nVar) >= 10

selection <- qsel & dpfilt 

print(paste("Filtered out: ", sum(!selection), sep=""))
print(paste("Variants selected: ", sum(selection), sep=""))
return(qc)
}

addQualTransCor <- function(qc, lociFilt) {


# lociFilt contains snp, cosmic and clinvar positions

qc$del <- as.character(as.vector(qc$del))
#qc$del <- as.character(as.vector(qc$del))
qc$ins <- as.character(as.vector(qc$ins))
qc$pVal <- 1
qc$pValAdj <- 1
#qc$pValPois <- 1
qc$qScore <- 0
qc$qScoreI <- 0
qc$qScoreD <- 0

pVal.I <- rep(1, nrow(qc))
pVal.D <- rep(1, nrow(qc))

rle(as.vector(qc$ntRef)) -> srle
hompols <- srle$lengths >= 5
homPosS <- cumsum(srle$lengths)[hompols] 
homPosE <- cumsum(srle$lengths)[which(hompols) - 1] 
homFilt <- rep(F, nrow(qc))
homFilt[homPosS] <- T
homFilt[homPosE] <- T

delSel <- qc$del != "." & !homFilt
if (sum(delSel > 0)) {
	

data.frame(matrix(unlist(strsplit(unlist(strsplit(qc$del[delSel], "\\|")), ":")), ncol=2, byrow=T)) -> dels
dels$ln <- nchar(as.vector(dels$X1))
colnames(dels) <- c("seq", "occ", "ln")
delPees <- ppois(1:50, mean(rep(dels$ln, dels$occ)), lower.tail=F)
}

insSel <- qc$ins != "." & !homFilt
if (sum(insSel) > 0) {
data.frame(matrix(unlist(strsplit(unlist(strsplit(qc$ins[insSel], "\\|")), ":")), ncol=2, byrow=T)) -> inss
inss$ln <- nchar(as.vector(inss$X1))
colnames(inss) <- c("seq", "occ", "ln")
insPees <- ppois(1:50, mean(rep(inss$ln, inss$occ)), lower.tail=F)
}
# Calculate transition bias


#calcTransCor(qc) -> pCor
globalProb <- sum(qc$nVar) / (sum(qc$nRef) + sum(qc$nVar))
# Calculate errorrate per amplicon
amps <- unique(qc$amp)
for (i in amps) {
	#print(i)
	sel <- qc$amp == i
        prob1 <- sum(qc$nVar[sel]) / (sum(qc$nRef[sel]) + sum(qc$nVar[sel])) 
       	prob <- prob1 * calcTransCor(qc[sel,], lociFilt)
	pees <- pbinom(qc$nVar[sel], qc$nRef[sel] + qc$nVar[sel], prob , lower.tail=F)
        qc$pVal[sel] <- pees
#       pois <- ppois(qc$nVar[sel], (qc$nRef + qc$nVar) * prob, lower.tail=F)
#       qc$pValPois[sel] <- pois

	# Seperate for indel/deletion events
	selDels <-  delSel & sel 
	if (sum(selDels) > 0) {
		data.frame(matrix(unlist(strsplit(unlist(strsplit(as.vector(qc$del[selDels]), "\\|")), ":")), ncol=2, byrow=T)) -> dels2
		dels2$ln <- nchar(as.vector(dels2$X1))
		colnames(dels2) <- c("seq", "occ", "ln")
		dels2$occ <- as.numeric(as.vector(dels2$occ))

	for (j in which(selDels)) {
		data.frame(matrix(unlist(strsplit(unlist(strsplit(as.vector(qc$del[j]), "\\|")), ":")), ncol=2, byrow=T)) -> dels
		dels$ln <- nchar(as.vector(dels$X1))
		colnames(dels) <- c("seq", "occ", "ln")
		dels$occ <- as.numeric(as.vector(dels$occ))
		#dels$pDel <- delPees[as.integer(dels$ln)]
		#dels$pVal <- ppois(as.numeric(dels$occ), dels$pDel, lower.tail=F)
		dels$pDel <- sapply(dels$ln, function(x) {sum(dels2$occ[dels2$ln == x]) / (sum(qc$nRef[sel]) + sum(qc$nVar[sel])) })
		dels$pVal <- pbinom(dels$occ, qc$dp[j], dels$pDel, lower.tail=F)
		dels$Q <- -10 * log10(dels$pVal)
			
		high <- order(dels$occ, -dels$Q, decreasing=T)[1]
		qc$del[j] <- paste(dels$seq[high], dels$occ[high], sep=":")		
		qc$qScoreD[j] <- dels$Q[high]
		pVal.D[j] <- dels$pVal[high]
	}
	}
	selIns <-  insSel & sel
	if (sum(selIns) > 0) {
		data.frame(matrix(unlist(strsplit(unlist(strsplit(as.vector(qc$ins[selIns]), "\\|")), ":")), ncol=2, byrow=T)) -> inss2
		inss2$ln <- nchar(as.vector(inss2$X1))
		colnames(inss2) <- c("seq", "occ", "ln")
		inss2$occ <- as.numeric(as.vector(inss2$occ))
	
	for (j in which(selIns)) {
		data.frame(matrix(unlist(strsplit(unlist(strsplit(as.vector(qc$ins[j]), "\\|")), ":")), ncol=2, byrow=T)) -> inss
		inss$ln <- nchar(as.vector(inss$X1))
		colnames(inss) <- c("seq", "occ", "ln")
		inss$occ <- as.numeric(as.vector(inss$occ))
		#inss$pIns <- insPees[as.integer(inss$ln)]
		#inss$pVal <- ppois(as.numeric(inss$occ), inss$pIns, lower.tail=F)
		inss$pIns <- sapply(inss$ln, function(x) {sum(inss2$occ[inss2$ln == x]) / (sum(qc$nRef[sel]) + sum(qc$nVar[sel])) })
		inss$pVal <- pbinom(inss$occ, qc$dp[j], inss$pIns, lower.tail=F)
		inss$Q <- -10 * log10(inss$pVal)

		high <- order(inss$occ, -inss$Q, decreasing=T)[1]
		qc$ins[j] <- paste(inss$seq[high], inss$occ[high], sep=":")
		qc$qScoreI[j] <- inss$Q[high]
		pVal.I[j] <- inss$pVal[high]
}	}
	#print(sum(-10 * log10(pees) > 100))
	#print(sum(qc$qScoreI[selIns] > 100))
	#print(sum(qc$qScoreD[selDels] > 100))
}

p.adjust(qc$pVal) -> qc$pValAdj
-10 * log(p.adjust(pVal.I)) -> qc$qScoreI
-10 * log(p.adjust(pVal.D)) -> qc$qScoreD
qc$qScore <- -10 * log10(qc$pValAdj)
#qc$qScore[qc$pValAdj == 0] <- 1000
qc$qScore[qc$qScore > 1000] <- 1000
qc$qScoreI[qc$qScoreI > 1000] <- 1000
qc$qScoreD[qc$qScoreD > 1000] <- 1000
qsel <- (qc$qScore >= 100 & qc$nVar >= 10) | (qc$qScoreI >= 200) | (qc$qScoreD >= 200)
dpfilt <- as.numeric(qc$dp) >= 100
#varfilt <- as.numeric(qc$nVar) >= 10

selection <- qsel & dpfilt 

print(paste("Filtered out: ", sum(!selection), sep=""))
print(paste("Variants selected: ", sum(selection), sep=""))
return(qc)
}


addQualVarCor <- function(qc) {
# Corrected for variant occurence

qc$del <- as.character(as.vector(qc$del))
qc$ins <- as.vector(qc$ins)
qc$pVal <- 1
qc$pValAdj <- 1
#qc$pValPois <- 1
qc$qScore <- 0
qc$qScoreI <- 0
qc$qScoreD <- 0

pVal.I <- rep(1, nrow(qc))
pVal.D <- rep(1, nrow(qc))

delSel <- qc$del != "."
data.frame(matrix(unlist(strsplit(unlist(strsplit(as.vector(qc$del[delSel]), "\\|")), ":")), ncol=2, byrow=T)) -> dels
dels$ln <- nchar(as.vector(dels$X1))
colnames(dels) <- c("seq", "occ", "ln")
delPees <- ppois(1:50, mean(rep(dels$ln, dels$occ)), lower.tail=F)

insSel <- qc$ins != "."
data.frame(matrix(unlist(strsplit(unlist(strsplit(as.vector(qc$ins[insSel]), "\\|")), ":")), ncol=2, byrow=T)) -> inss
inss$ln <- nchar(as.vector(inss$X1))
colnames(inss) <- c("seq", "occ", "ln")
insPees <- ppois(1:50, mean(rep(dels$ln, dels$occ)), lower.tail=F)

# Calculate base bias per run
sum(qc$nA[grepl("A", qc$ntRef, ignore.case=T)]) -> refA
sum(qc$nA[!grepl("A", qc$ntRef, ignore.case=T)]) -> varA
sum(qc$nC[grepl("C", qc$ntRef, ignore.case=T)]) -> refC
sum(qc$nC[!grepl("C", qc$ntRef, ignore.case=T)]) ->varC
sum(qc$nG[grepl("G", qc$ntRef, ignore.case=T)]) -> refG
sum(qc$nG[!grepl("G", qc$ntRef, ignore.case=T)]) -> varG
sum(qc$nT[grepl("T", qc$ntRef, ignore.case=T)]) -> refT
sum(qc$nT[!grepl("T", qc$ntRef, ignore.case=T)]) -> varT

varA / (refA + varA) -> errA
varC / (refC + varC) -> errC
varG / (refG + varG) -> errG
varT / (refT + varT) -> errT
c(errA, errC, errG, errT) / mean(c(errA, errC, errG, errT)) -> errACGT

errMat <- qc[,9:12]
errMat[grepl("A", qc$ntRef, ignore.case=T),1] <- 0
errMat[grepl("C", qc$ntRef, ignore.case=T),2] <- 0
errMat[grepl("G", qc$ntRef, ignore.case=T),3] <- 0
errMat[grepl("T", qc$ntRef, ignore.case=T),4] <- 0

pCor <- rep(1, nrow(qc))
pCor[grepl("^A", qc$ntVar, ignore.case=T)] <- errACGT[1]
pCor[grepl("^C", qc$ntVar, ignore.case=T)] <- errACGT[2]
pCor[grepl("^G", qc$ntVar, ignore.case=T)] <- errACGT[3]
pCor[grepl("^T", qc$ntVar, ignore.case=T)] <- errACGT[4]

# Calculate errorrate per amplicon
amps <- unique(qc$amp)
for (i in amps) {
	#print(i)
	sel <- qc$amp == i
        prob <- sum(qc$nVar[sel]) / (sum(qc$nRef[sel]) + sum(qc$nVar[sel]))
	#prob <- median(qc$nVar[sel] / (qc$nRef[sel] + qc$nVar[sel]))
        pees <- pbinom(qc$nVar[sel], qc$nRef[sel] + qc$nVar[sel], prob * pCor[sel], lower.tail=F)
        qc$pVal[sel] <- pees
#       pois <- ppois(qc$nVar[sel], (qc$nRef + qc$nVar) * prob, lower.tail=F)
#       qc$pValPois[sel] <- pois

	# Seperate for indel/deletion events
	selDels <-  delSel & sel
	if (sum(selDels) > 0) {
		data.frame(matrix(unlist(strsplit(unlist(strsplit(as.vector(qc$del[selDels]), "\\|")), ":")), ncol=2, byrow=T)) -> dels2
		dels2$ln <- nchar(as.vector(dels2$X1))
		colnames(dels2) <- c("seq", "occ", "ln")
		dels2$occ <- as.numeric(as.vector(dels2$occ))
	}
	for (j in which(selDels)) {
		data.frame(matrix(unlist(strsplit(unlist(strsplit(as.vector(qc$del[j]), "\\|")), ":")), ncol=2, byrow=T)) -> dels
		dels$ln <- nchar(as.vector(dels$X1))
		colnames(dels) <- c("seq", "occ", "ln")
		dels$occ <- as.numeric(as.vector(dels$occ))
		#dels$pDel <- delPees[as.integer(dels$ln)]
		#dels$pVal <- ppois(as.numeric(dels$occ), dels$pDel, lower.tail=F)
		dels$pDel <- sapply(dels$ln, function(x) {sum(dels2$occ[dels2$ln == x]) / (sum(qc$nRef[sel]) + sum(qc$nVar[sel])) })
		dels$pVal <- pbinom(dels$occ, qc$dp[j], dels$pDel, lower.tail=F)
		dels$Q <- -10 * log10(dels$pVal)
			
		high <- order(dels$occ, -dels$Q, decreasing=T)[1]
		qc$del[j] <- paste(dels$seq[high], dels$occ[high], sep=":")		
		qc$qScoreD[j] <- dels$Q[high]
		pVal.D[j] <- dels$pVal[high]
	}
	selIns <-  insSel & sel
	if (sum(selIns) > 0) {
		data.frame(matrix(unlist(strsplit(unlist(strsplit(as.vector(qc$ins[selIns]), "\\|")), ":")), ncol=2, byrow=T)) -> inss2
		inss2$ln <- nchar(as.vector(inss2$X1))
		colnames(inss2) <- c("seq", "occ", "ln")
		inss2$occ <- as.numeric(as.vector(inss2$occ))
	}
	for (j in which(selIns)) {
		data.frame(matrix(unlist(strsplit(unlist(strsplit(as.vector(qc$ins[j]), "\\|")), ":")), ncol=2, byrow=T)) -> inss
		inss$ln <- nchar(as.vector(inss$X1))
		colnames(inss) <- c("seq", "occ", "ln")
		inss$occ <- as.numeric(as.vector(inss$occ))
		#inss$pIns <- insPees[as.integer(inss$ln)]
		#inss$pVal <- ppois(as.numeric(inss$occ), inss$pIns, lower.tail=F)
		inss$pIns <- sapply(inss$ln, function(x) {sum(inss2$occ[inss2$ln == x]) / (sum(qc$nRef[sel]) + sum(qc$nVar[sel])) })
		inss$pVal <- pbinom(inss$occ, qc$dp[j], inss$pIns, lower.tail=F)
		inss$Q <- -10 * log10(inss$pVal)

		high <- order(inss$occ, -inss$Q, decreasing=T)[1]
		qc$ins[j] <- paste(inss$seq[high], inss$occ[high], sep=":")
		qc$qScoreI[j] <- inss$Q[high]
		pVal.I[j] <- inss$pVal[high]
	}
	#print(sum(-10 * log10(pees) > 100))
	#print(sum(qc$qScoreI[selIns] > 100))
	#print(sum(qc$qScoreD[selDels] > 100))
}

p.adjust(qc$pVal) -> qc$pValAdj
-10 * log(p.adjust(pVal.I)) -> qc$qScoreI
-10 * log(p.adjust(pVal.D)) -> qc$qScoreD
qc$qScore <- -10 * log10(qc$pValAdj)
#qc$qScore[qc$pValAdj == 0] <- 1000
qc$qScore[qc$qScore > 1000] <- 1000
qc$qScoreI[qc$qScoreI > 1000] <- 1000
qc$qScoreD[qc$qScoreD > 1000] <- 1000
qsel <- (qc$qScore >= 100 & qc$nVar >= 10) | (qc$qScoreI >= 200) | (qc$qScoreD >= 200)
dpfilt <- as.numeric(qc$dp) >= 100
#varfilt <- as.numeric(qc$nVar) >= 10

selection <- qsel & dpfilt 

print(paste("Filtered out: ", sum(!selection), sep=""))
print(paste("Variants selected: ", sum(selection), sep=""))
return(qc)
}

#addQual(qc) -> qc
#write.table(qc, output, row.names=F, col.names=T, quote=F, sep="\t")
#write.table(qc[selection,], output, row.names=F, col.names=T, quote=F, sep="\t")

calcTrans <- function(qcRaw, lociFilt = NULL) {

sel <- rep(T, nrow(qcRaw))
if (!is.null(lociFilt)) {
paste(qcRaw$X.chr, qcRaw$pos, sep=":") -> qcpos
sel <- ! qcpos %in% lociFilt
}

qc <- qcRaw[sel,]

# Calculate transition bias
refA <- grepl("A", qc$ntRef, ignore.case=T)
refC <- grepl("C", qc$ntRef, ignore.case=T)
refG <- grepl("G", qc$ntRef, ignore.case=T)
refT <- grepl("T", qc$ntRef, ignore.case=T)
varA <- grepl("A", qc$ntVar, ignore.case=T)
varC <- grepl("C", qc$ntVar, ignore.case=T)
varG <- grepl("G", qc$ntVar, ignore.case=T)
varT <- grepl("T", qc$ntVar, ignore.case=T)
varN <- grepl(".", qc$ntVar, ignore.case=T)
callA <- grepl("^A", qc$ntVar, ignore.case=T)
callC <- grepl("^C", qc$ntVar, ignore.case=T)
callG <- grepl("^G", qc$ntVar, ignore.case=T)
callT <- grepl("^T", qc$ntVar, ignore.case=T)
ntA <- sum(refA)
ntC <- sum(refC)
ntG <- sum(refG)
ntT <- sum(refT)

callList <- list(callA, callC, callG, callT)
varList <- list(varA, varC, varG, varT)
refList <- list(refA, refC, refG, refT)

#tots <- c( 
#rep(sum(qc$nRef[refA], qc$nVar[refA]), 4),
#rep(sum(qc$nRef[refC], qc$nVar[refC]), 4),
#rep(sum(qc$nRef[refG], qc$nVar[refG]), 4),
#rep(sum(qc$nRef[refT], qc$nVar[refT]), 4))

nvar <- c( 
rep(sum(qc$nC[refA], qc$nG[refA], qc$nT[refA]), 4),
rep(sum(qc$nA[refC], qc$nG[refC], qc$nT[refC]), 4),
rep(sum(qc$nA[refG], qc$nC[refG], qc$nT[refG]), 4),
rep(sum(qc$nA[refT], qc$nC[refT], qc$nG[refT]), 4))

ncall <- c( 
rep(sum(qc$nVar[refA]), 4),
rep(sum(qc$nVar[refC]), 4),
rep(sum(qc$nVar[refG]), 4),
rep(sum(qc$nVar[refT]), 4))

nref <- c( 
rep(sum(qc$nRef[refA]), 4),
rep(sum(qc$nRef[refC]), 4),
rep(sum(qc$nRef[refG]), 4),
rep(sum(qc$nRef[refT]), 4))

tots <- nref + nvar

nts <- c("A", "C", "G", "T")
sumTrans <- data.frame(trans=paste(rep(nts, each=4), rep(nts, 4), sep=":"), npos=rep(c(ntA, ntC, ntG, ntT), each=4), cnt=rep(0, 16), tot=tots, nvar=nvar, nref=nref, ncall=ncall)

for (i in 1:4) {
        pos <- (i - 1) * 4
for (j in 1:4) {
#       print(sum(refList[[i]] & varList[[j]]))

        cnt <- sum(qc[refList[[i]] & varList[[j]], 8 + j])
        if (i == j) {
#               cnt <- sum(qc[refList[[i]] & varN, 8 + j])
                cnt <- NA
        }
        sumTrans$cnt[pos + j] <- cnt
        #print(cnt)
}}

#sumTrans$cnt / sumTrans$nvar / sumTrans$tot * 1000000 -> tmp
sumTrans$cnt / sumTrans$tot * 1000000 -> tmp

#totalNT <- sum(qc$nVar) + sum(qc$nRef)
#sumTrans$cnt / median(sumTrans$cnt, na.rm=T) -> sumTrans$cor
sumTrans$cnt / median(sumTrans$cnt, na.rm=T) -> sumTrans$cor
tmp / median(tmp, na.rm=T) -> sumTrans$cor2

sumTrans$cor[sumTrans$npos <= 15] <- 1
sumTrans$cor2[sumTrans$npos <= 15] <- 1

#sumTrans$cnt / totalNT -> sumTrans$cor
#sumTrans$cnt / min(sumTrans$cnt, na.rm=T) -> sumTrans$cor

return(sumTrans)

}

calcTransCor <- function(qcRaw, lociFilt = NULL) {

calcTrans(qcRaw, lociFilt) -> sumTrans

pCor <- rep(1, nrow(qcRaw))
paste(qcRaw$ntRef, substr(qcRaw$ntVar, 1, 1), sep=":") -> qcTrans

for (i in c("A", "C", "G", "T")) {
for (j in c("A", "C", "G", "T")) {
        if (i != j) {
		tTrans <- paste(i, j, sep=":")
		pCor[grepl(tTrans, qcTrans, ignore.case=T)] <- sumTrans$cor2[sumTrans$trans == tTrans]
        }

}}
return(pCor)
}


calcTransCorOld <- function(qc) {
# Calculate transition bias
refA <- grepl("A", qc$ntRef, ignore.case=T)
refC <- grepl("C", qc$ntRef, ignore.case=T)
refG <- grepl("G", qc$ntRef, ignore.case=T)
refT <- grepl("T", qc$ntRef, ignore.case=T)
varA <- grepl("A", qc$ntVar, ignore.case=T)
varC <- grepl("C", qc$ntVar, ignore.case=T)
varG <- grepl("G", qc$ntVar, ignore.case=T)
varT <- grepl("T", qc$ntVar, ignore.case=T)
varN <- grepl(".", qc$ntVar, ignore.case=T)
callA <- grepl("^A", qc$ntVar, ignore.case=T)
callC <- grepl("^C", qc$ntVar, ignore.case=T)
callG <- grepl("^G", qc$ntVar, ignore.case=T)
callT <- grepl("^T", qc$ntVar, ignore.case=T)

callList <- list(callA, callC, callG, callT)
varList <- list(varA, varC, varG, varT)
refList <- list(refA, refC, refG, refT)

tots <- c( 
rep(sum(qc$nRef[refA], qc$nVar[refA]), 4),
rep(sum(qc$nRef[refC], qc$nVar[refC]), 4),
rep(sum(qc$nRef[refG], qc$nVar[refG]), 4),
rep(sum(qc$nRef[refT], qc$nVar[refT]), 4))

nvar <- c( 
rep(sum(qc$nVar[refA]), 4),
rep(sum(qc$nVar[refC]), 4),
rep(sum(qc$nVar[refG]), 4),
rep(sum(qc$nVar[refT]), 4))

nref <- c( 
rep(sum(qc$nRef[refA]), 4),
rep(sum(qc$nRef[refC]), 4),
rep(sum(qc$nRef[refG]), 4),
rep(sum(qc$nRef[refT]), 4))

nts <- c("A", "C", "G", "T")
sumTrans <- data.frame(trans=paste(rep(nts, each=4), rep(nts, 4), sep=":"), cnt=rep(0, 16), tot=tots, nvar=nvar, nref=nref)

for (i in 1:4) {
        pos <- (i - 1) * 4
for (j in 1:4) {
#       print(sum(refList[[i]] & varList[[j]]))

        cnt <- sum(qc[refList[[i]] & varList[[j]], 8 + j])
        if (i == j) {
#               cnt <- sum(qc[refList[[i]] & varN, 8 + j])
                cnt <- NA
        }
        sumTrans$cnt[pos + j] <- cnt
        #print(cnt)
}}
#totalNT <- sum(qc$nVar) + sum(qc$nRef)
sumTrans$cnt / median(sumTrans$cnt, na.rm=T) -> sumTrans$cor
#sumTrans$cnt / totalNT -> sumTrans$cor
#sumTrans$cnt / min(sumTrans$cnt, na.rm=T) -> sumTrans$cor
#print(sumTrans)

pCor <- rep(1, nrow(qc))
for (i in 1:4) {
        pos <- (i - 1) * 4
for (j in 1:4) {
        if (i != j) {
                pCor[refList[[i]] & callList[[j]]] <- sumTrans$cor[pos + j]
        }

}}
return(pCor)
}
