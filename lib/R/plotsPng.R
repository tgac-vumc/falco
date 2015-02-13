#!/usr/bin/Rscript 

args <- commandArgs(T)

print(args[1]) # qc.ann.qual.txt
print(args[2]) # qc2.ann.qtxt
print(args[3]) # qc.targets.txt
print(args[4]) # res.txt (unfiltered)
print(args[5]) # clinvar
print(args[6]) # lociFilt
print(args[7]) # basename
options(warn=1)
print(paste("Reading", args[1], sep=" "));
read.delim(args[1], header=T, stringsAsFactors=F) -> d
print(paste("Reading", args[2], sep=" "));
read.delim(args[2], header=T, stringsAsFactors=F) -> d2
print(paste("Reading", args[3], sep=" "));
read.delim(args[3], header=T, stringsAsFactors=F) -> qc
print(paste("Reading", args[4], sep=" "));
read.delim(args[4], header=T, stringsAsFactors=F, row.names=NULL) -> res
print(paste("Reading", args[5], sep=" "));
read.delim(args[5], header=T, stringsAsFactors=F) -> dclin
print(paste("Reading", args[6], sep=" "));
read.delim(args[6], header=F, stringsAsFactors=F) -> lociFilt
base <- args[7]

# Error rate per cycle
# VAF distribution  
# varDepth distribution ?
# Print coverage plots per amplicon

tmpCP <- paste(d$X.chr, d$pos, sep=":")
tmpClin <- paste(paste("chr", dclin[,1], sep=""), dclin[,2], sep=":")
tmpSnp <- paste(paste("chr", lociFilt[,1], sep=""), lociFilt[,2], sep=":")

inClin <- tmpCP %in% tmpClin 

inSnp <- tmpCP %in% tmpSnp

dir.create(base)
setwd(base)
#unique(unlist(strsplit(d$V3, "\\|"))) -> targets
targets <- qc[,1]
if (1) { 
print("Plotting variant heatmap")
# Including off targets increases the amplicon length set at median amplicon length + 10%.
median(qc$end - qc$start) -> maxAmpLn
maxAmpLn <- ceiling(1.2 * maxAmpLn)
print(maxAmpLn)
ampHeatF <- matrix(ncol=maxAmpLn, nrow=nrow(qc))
ampHeatR <- matrix(ncol=maxAmpLn, nrow=nrow(qc))
#dPos <- paste(d$chr, d$pos, sep=":")
#errPos <- paste(dclin$chr, dclin$pos, sep=":")
#dPosClin <- dPos %in% errPos

dPosClin <- !inClin & !inSnp 

for (i in 1:nrow(qc)) {
	d$pos %in% qc$start[i]:qc$end[i] -> sel
	poss <- d$pos[sel & !dPosClin] - qc$start[i] + 1
	poss <- poss[poss <= maxAmpLn]
	ampHeatF[i, poss] <- d$nVar[sel & !dPosClin] + d$nN[sel & !dPosClin]
	ampHeatR[i, poss] <- rev(d$nVar[sel & !dPosClin] + d$nN[sel & !dPosClin])
}
png(paste(base,"heat","png", sep="."), width=960, height=480)
par(mar=c(4,4,4,2) + .1)

ampHeatF[is.na(ampHeatF)] <- 0
ampHeatR[is.na(ampHeatR)] <- 0
layout(matrix(1:2, nrow=2))
maxN <- quantile(ampHeatF, .99)
#heatmap(ampHeat)
boxplot(ampHeatF, ylim=c(0, maxN), pch=20, cex=.4, main="Non reference counts R1")
rect(150, 0, ncol(ampHeatF), maxN, col="#00000040")
boxplot(ampHeatR, ylim=c(0, maxN), pch=20, cex=.4, main="Non reference counts R2")
rect(150, 0, ncol(ampHeatR), maxN, col="#00000040")
dev.off()

#d[,7]/d[,4] * 100 -> d$pctREF
#d[,8]/d[,4] * 100 -> d$pctVAR
d$nRef/d$dp * 100 -> d$pctREF
d$nVar/d$dp * 100 -> d$pctVAR
d$vaf <- d$nVar / (d$nVar + d$nRef) * 100
print("Plotting ref frequencies")
png(paste(base,"raf","png", sep="."), width=480, height=480)
dref <- density(d$pctREF[!dPosClin])
#plot(dref, ylim=range(dref$y[dref$x < 80] ))
hist(d$pctREF[!dPosClin], breaks=100, ylim=c(0, 30), col="#00000040", border=NA)
dev.off()

print("Plotting var frequencies")
png(paste(base,"vaf","png", sep="."), width=480, height=480)
dvar <- density(d$pctVAR[!dPosClin])
#plot(dvar, ylim=range(dvar$y[dvar$x > 20] ))
hist(d$pctVAR[!dPosClin], breaks=100, ylim=c(0, 30), col="#00000040", border=NA)
dev.off()

# Set vaf cutoff
#q95 <- quantile(d$qScore, .95, na.rm=T)
q95 <- 20
vs95 <- sd(d$pctVAR[!dPosClin & d$qScore < q95], na.rm=T)
v95 <- mean(d$pctVAR[!dPosClin & d$qScore < q95], na.rm=T) + vs95 * 3

nCall <- sum((d$pctVAR > v95) & d$qScore > q95 & !inSnp, na.rm=T)
print(nCall)
print("Plotting vaf vs q")
png(paste(base,"snv-q","png", sep="."), width=480, height=480)
plot(d$qScore, d$pctVAR, pch=20, cex=.5, main=paste("q:", q95, "-", "v:", v95, "--", "s:", nCall, sep=" "))
abline(h=c(1,5), v=c(20,13), lty=2)
abline(v=q95, h=v95, lty=3, col=2)
points(d$qScore[d$dp < 100], d$pctVAR[d$dp < 100], pch=20, cex=1.2, col=2)
points(d$qScore[d$nVar < 10], d$pctVAR[d$nVar < 10], pch=20, cex=1, col=4)
legend(30,80, c("pass", "dp < 100", "nVar < 10"), col=c(1,2,4), pch=20)
dev.off()

print("Plotting vaf vs q zoomed")
png(paste(base,"snv-q-zoom","png", sep="."), width=480, height=480)
plot(d$qScore, d$pctVAR, pch=20, cex=.5, xlim=c(0, 40),ylim=c(0,10))
abline(h=c(1,5), v=c(20,13), lty=2)
abline(v=q95, h=v95, lty=3, col=2)
points(d$qScore[d$dp < 100], d$pctVAR[d$dp < 100], pch=20, cex=1.2, col=2)
points(d$qScore[d$nVar < 10], d$pctVAR[d$nVar < 10], pch=20, cex=1, col=4)
legend("topright", c("pass", "dp < 100", "nVar < 10"), col=c(1,2,4), pch=20)
dev.off()

print("Plotting ins vs q")
d$vafI <- 0
d$occI <- 0
if (sum(d$ins != ".") > 0) {
png(paste(base,"ins-q","png", sep="."), width=480, height=480)
data.frame(matrix(unlist(strsplit(unlist(strsplit(as.character(d$ins[d$ins != "."]), "\\|")), ":")), ncol=2, byrow=T)) -> insDat
colnames(insDat) <- c("seq", "occ")
insDat$seq <- as.character(insDat$seq)
insDat$occ <- as.numeric(as.vector(insDat$occ))
nchar(insDat$seq) -> insDat$ln
d$occI[d$ins != "."] <- insDat$occ
d$vafI[d$ins != "."] <- insDat$occ / d$nRef[d$ins != "."]
d$qScoreI[d$qScoreI > 1000] <- 1000
plot(d$qScoreI, d$occI, pch=20, cex=1)
abline(h=c(10), v=c(20,13), lty=2)
# dp 500
sel <- d$dp < 100
points(d$qScoreI[sel], d$occI[sel], pch=20, cex=1.2, col=2)
# nVar 10
#sel <- d$occI < 10
#points(d$qScoreI[sel], d$occI[sel], pch=20, cex=1, col=3)
# vaf < 0.01
sel <- d$vafI < 0.01
points(d$qScoreI[sel], d$occI[sel], pch=20, cex=1, col=4)
legend("topleft", c("pass", "dp < 100", "vaf < 0.01"), col=c(1,2,4), pch=20)
dev.off()
}

print("Plotting del vs q")
d$vafD <- 0
d$occD <- 0
if (sum(d$del != "." ) > 0) {
png(paste(base,"del-q","png", sep="."), width=480, height=480)
data.frame(matrix(unlist(strsplit(unlist(strsplit(as.character(d$del[d$del != "."]), "\\|")), ":")), ncol=2, byrow=T)) -> delDat
colnames(delDat) <- c("seq", "occ")
delDat$seq <- as.character(delDat$seq)
delDat$occ <- as.numeric(as.vector(delDat$occ))
nchar(delDat$seq) -> delDat$ln
d$occD[d$del != "."] <- delDat$occ
d$vafD[d$del != "."] <- d$occD[d$del != "."] / d$nRef[d$del != "."]
d$qScoreD[d$qScoreD > 1000] <- 1000
#hist(d$qScoreD, breaks=1000, ylim=c(0,100))
#print(sum(d$qScoreD > 20))
#Sys.sleep(1)
plot(d$qScoreD, d$occD, pch=20, cex=1)
abline(h=c(10), v=c(20,13), lty=2)
# dp 500
sel <- d$dp < 100
points(d$qScoreD[sel], d$occD[sel], pch=20, cex=1.2, col=2)
# nVar 10
#sel <- d$occD < 10
#points(d$qScoreD[sel], d$occD[sel], pch=20, cex=1, col=3)
# vaf < 0.01
sel <- d$vafD < 0.01
points(d$qScoreD[sel], d$occD[sel], pch=20, cex=1, col=4)
legend("topleft", c("pass", "dp < 100", "vaf < 0.01"), col=c(1,2,4), pch=20)
dev.off()
}

sub("chr", "", qc$chr) -> qc$chrn
ampOrd <- order(qc$chrn, qc$start)

print("Plotting amplicon depths")
png(paste(base,"amp-dp","png", sep="."), width=960, height=480)
barplot(qc$depth[ampOrd])
abline(h=100, lty=2)
totRead <- sum(qc$depth, na.rm=T)
legend("topright", paste("Total read count: ", totRead, sep=""))
dev.off()
}
#
res$POS <- as.numeric(as.vector(res$POS))
d$pos <- as.numeric(as.vector(d$pos))

grepl("a", d$ntRef, ignore.case=T) -> Asel
grepl("c", d$ntRef, ignore.case=T) -> Csel
grepl("g", d$ntRef, ignore.case=T) -> Gsel
grepl("t", d$ntRef, ignore.case=T) -> Tsel
cntr <- 0

biasTable <- data.frame(matrix(NA, nrow=length(targets),  ncol=12), row.names=targets)
colnames(biasTable) <- c("AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG")

for (i in targets) {
	sel <- grepl(i, d$amp)
	print(i)
	for (c in unique(d[sel,1])) {
		csel <- d[,1] == c
		print(c)
		reg <- paste(range(d[sel & csel,2]), collapse="-", sep="-")
		creg <- paste(c, reg, sep=":")
	
		qcsel <- qc[,1] == i	
		allsel <- sel & csel
		clinSel <- dclin$amp == i
		clinPos <- dclin$pos[clinSel]
#		inClin <- d$pos %in% clinPos	
	
	        #dels <- (d[,15] != ".") & allsel
	        #ins <- (d[,14] != ".") & sel	
	#	print(reg)
	#	print(creg)

		# Select annotations
		# calculate quantiles and plot them
#		quantile(d$pctVAR[sel & csel], c(.75,.90,.95,.99), na.rm=T) -> cutoffs
		tmp <- d[sel & csel,]
		#dels <- tmp[,15] != "."
		#ins <- tmp[,14] != "."
		dels <- allsel & d$del != "."
		ins <- allsel & d$ins != "."

		resTmp <- res[(res$POS >= min(tmp[,2])) & (res$POS <= max(tmp[,2])) & (res$X.CHROM == c),]
		#resSel <- res$TARGET == i	
		png(paste(base, i, "cov", "png", sep="."), width=960, height=480)
		
		xvec <- c(0, rep(1:(sum(allsel) - 1), each=2), sum(allsel))
		layout(matrix(1:2, ncol=1), heights=c(2,1))
		par(mar=c(0, 4, 4, 7) + .1, xaxt="n")
		print("Cov")	
		plot(NULL, type="n", main=i, xlab="Position", ylab="Depth", xlim=range(d$pos[allsel]), ylim=range(0, d$dp[allsel]))
		xvec <- d$pos[allsel]
		yvec <- d$dp[allsel]
	
		polygon(rep(c(xvec, max(xvec) + 1), each=2) - .5, c(0, rep(yvec, each=2), 0), col="#00000010", border=1)
		polygon(rep(c(xvec, max(xvec) + 1), each=2) - .5, c(0, rep(rowSums(d[allsel, 7:8]), each=2), 0), col="#00000010")
#		polygon(rep(c(xvec, max(xvec) + 1), each=2) - .5, c(0, rep(d$nVar[allsel], each=2), 0), col=2)
		varCol <- rep(2, sum(allsel))
		varCol[grepl("^a", d$ntVar[allsel], ignore.case=T)] <- "green"
		varCol[grepl("^c", d$ntVar[allsel], ignore.case=T)] <- "blue"
		varCol[grepl("^g", d$ntVar[allsel], ignore.case=T)] <- "black"
		varCol[grepl("^t", d$ntVar[allsel], ignore.case=T)] <- "red"
		rect(xvec-.5,d$nVar[allsel],xvec+.5,0, col=varCol)


		ampMax <- max(d[allsel,4])
		#abline(v=resTmp$POS, lty=2, col=2)
		#delMark <- list()
		#insMark <- list()
			
		for (idx in which(dels)) {
			matrix(unlist(strsplit(as.character(as.vector(d[idx,15])), "[:|]")), ncol=2, byrow=T) -> delDat
				nchar(delDat[,1]) -> lngts
				dpths <- as.numeric(delDat[,2])
				x1 <- d[idx,2] + .5 
				x2 <- x1 + lngts - 1
				if (x2 > max(d$pos[allsel])) {
					next
				}
		#	delMark[[length(delMark) + 1]] <- c(x1, x2)
			rect(x1, ampMax - dpths , x2, ampMax, col="#00000040")
		}

		for (idx in which(ins)) {
			matrix(unlist(strsplit(as.character(as.vector(d[idx,14])), "[:|]")), ncol=2, byrow=T) -> insDat
				nchar(insDat[,1]) -> lngts
				dpths <- as.numeric(insDat[,2])
				x1 <- d$pos[idx] + 1
		#		insMark[length(insMark) + 1] <- x1
				for (ii in 1:ncol(insDat)) {
					polygon(
							c(x1, x1 - (.5 * lngts[ii]), x1 + (.5 * lngts[ii])),
							c(ampMax, ampMax - dpths[ii], ampMax - dpths[ii]),
							, col="#00000040")
				}
		}
		par(xpd=NA)	
		legend(max(d$pos[allsel]) + sum(allsel) * 0.05, ampMax, legend=c("Filtered nt", "Ref nt", "Non-ref nt"), col=c("#00000040", "#00000080", "red"), pch=15)
		legend(max(d$pos[allsel]) + sum(allsel) * 0.05, 0, legend=c("A", "C", "G", "T"), col=c("green", "blue", "black", "red"), pch=15)
		par(xpd=F)	
		par(mar=c(5, 4, 0, 7) + .1, xaxt="s")
		plot(NULL, type="n", xlab="Position", ylab=NA, xlim=range(d$pos[allsel]), ylim=c(0,4.5), yaxt="n")
		axis(side=2, at=1:4, labels=c("ins", "del", "snv", "ref"), las=2)
		abline(h=1:4, col="#00000040")
		sres <- nchar(as.character(resTmp$REF)) == nchar(as.character(resTmp$ALT))
		dres <- nchar(resTmp$REF) > nchar(resTmp$ALT)
		ires <- nchar(resTmp$REF) < nchar(resTmp$ALT)
		
		resTmp$marks <- rep(3.4, nrow(resTmp))
		
		resTmp$cols <- rep("#00000040", nrow(resTmp))
		resTmp$cols[grepl("^a", resTmp$ALT, ignore.case=T) & sres] <- "green"
		resTmp$cols[grepl("^c", resTmp$ALT, ignore.case=T) & sres] <- "blue"
		resTmp$cols[grepl("^g", resTmp$ALT, ignore.case=T) & sres] <- "black"
		resTmp$cols[grepl("^t", resTmp$ALT, ignore.case=T) & sres] <- "red"

		#resTmp$cols[dres] <- 2
		#resTmp$cols[ires] <- 2
		resTmp$marks[dres] <- 2.4
		resTmp$marks[ires] <- 1.4
		resTmp$markP1 <- resTmp$POS -.5
		resTmp$markP2 <- resTmp$POS +.5
		resTmp$markP1[dres] <- resTmp$POS[dres] +.5
		resTmp$markP2[dres] <- resTmp$POS[dres] +.5 + nchar(resTmp$REF[dres]) - 1
		resTmp$markP1[ires] <- resTmp$POS[ires] +.5
		resTmp$markP2[ires] <- resTmp$POS[ires] +.5 + 1
		ntCol <- rep(1, sum(allsel))
		ntCol[grepl("a", d$ntRef[allsel], ignore.case=T)] <- "green"
		ntCol[grepl("c", d$ntRef[allsel], ignore.case=T)] <- "blue"
		ntCol[grepl("g", d$ntRef[allsel], ignore.case=T)] <- "black"
		ntCol[grepl("t", d$ntRef[allsel], ignore.case=T)] <- "red"
	
		rect(d$pos[allsel] -.5, 4.4, d$pos[allsel] +.5, 3.6, col=ntCol, border=0)
		rect(resTmp$markP1, resTmp$marks, resTmp$markP2, resTmp$marks - .8, col=resTmp$cols, border=0)	
		
		if(length(clinPos != 0)) {
			rect(clinPos - .5, 3.5, clinPos + .5, .5, border="black", lwd=2)
		}
		#lapply(delMark, function(x) { rect(x[1], 2, x[2], 1, col=2, border=NA)  })
		#lapply(insMark, function(x) { rect(x[1], 1, x[1] + 1, 0, col=3, border=NA)  })
		
		# Highlight clinically relevant loci	
		#for (i in which(dclin$amp == i)) {
		#	if (dclin$REF[i] == nchar
		#}
		
		dev.off()

		sansClin <- allsel & !inClin
		sansSnp <- allsel & !inSnp

		d$nA[sansClin & sansSnp & Asel] -> AA
		d$nC[sansClin & sansSnp & Asel] -> AC
		d$nG[sansClin & sansSnp & Asel] -> AG
		d$nT[sansClin & sansSnp & Asel] -> AT

		d$nA[sansClin & sansSnp & Csel] -> CA
		d$nC[sansClin & sansSnp & Csel] -> CC
		d$nG[sansClin & sansSnp & Csel] -> CG
		d$nT[sansClin & sansSnp & Csel] -> CT

		d$nA[sansClin & sansSnp & Gsel] -> GA
		d$nC[sansClin & sansSnp & Gsel] -> GC
		d$nG[sansClin & sansSnp & Gsel] -> GG
		d$nT[sansClin & sansSnp & Gsel] -> GT

		d$nA[sansClin & sansSnp & Tsel] -> TA
		d$nC[sansClin & sansSnp & Tsel] -> TC
		d$nG[sansClin & sansSnp & Tsel] -> TG
		d$nT[sansClin & sansSnp & Tsel] -> TT

		print("Bias")	
		png(paste(base, i, "bias", "png", sep="."), width=960, height=480)

		Acnt <- sum(sum(AA), sum(AC), sum(AG), sum(AT))
		Ccnt <- sum(sum(CA), sum(CC), sum(CG), sum(CT))
		Gcnt <- sum(sum(GA), sum(GC), sum(GG), sum(GT))
		Tcnt <- sum(sum(TA), sum(TC), sum(TG), sum(TT))

		NTsums <- c(
			sum(AC)/Acnt, sum(AG)/Acnt, sum(AT)/Acnt, 
			sum(CA)/Ccnt, sum(CG)/Ccnt, sum(CT)/Ccnt,
			sum(GA)/Gcnt, sum(GC)/Gcnt, sum(GT)/Gcnt,
			sum(TA)/Tcnt, sum(TC)/Tcnt, sum(TG)/Tcnt)
		NTsums[!is.finite(NTsums)] <- 0
		NTsums[is.nan(NTsums)] <- 0 
		NTsums[is.na(NTsums)] <- 0 
	
		barplot(NTsums, names.arg=c("AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"), main=i)
		dev.off()
		biasTable[rownames(biasTable) == i, ] <- NTsums

		print("Cov2")
		png(paste(base, i, "cov2", "png", sep="."), width=960, height=480)
		layout(matrix(1:2, ncol=1), heights=c(2,1))
		par(mar=c(0, 4, 4, 7) + .1, xaxt="n")
		

		dat <- d2[d2$amp == i,]
		rled <- rle(sort(dat$pos))	
		rds <- sort(unique(dat$X.read))
		lst <- c()
		for (ii in order(rled$lengths, decreasing=T)) 
		{

				rled$values[ii] -> pos
				dat$X.read[dat$pos == pos]
				lst <- c(lst, rds %in% dat$X.read[dat$pos == pos])

		}
		#xlim <- range(dat$pos)
		xlim=range(d$pos[allsel])
		ylim=range(0, d$dp[allsel])
		plot(NULL, xlim=xlim, ylim=ylim, main=i, ylab="Depth")
		rect(xlim[1] - 0.5, 0, xlim[2] + .5, ylim[2], col="#00000010")
		if (length(lst) > 0) {
		matrix(lst, nrow=length(rds), byrow=F) -> mat
		do.call( order, data.frame(mat) ) -> ord
		for (e in 1:length(ord)) {
				sel <- dat$X.read == rds[ord[e]]
				colvec <- rep(1, sum(sel))
				colvec[grepl("A$", dat$event[sel], ignore.case=T)] <- "green"
				colvec[grepl("C$", dat$event[sel], ignore.case=T)] <- "blue"
				colvec[grepl("G$", dat$event[sel], ignore.case=T)] <- "black"
				colvec[grepl("T$", dat$event[sel], ignore.case=T)] <- "red"
				#points(dat$pos[ sel ], rep(e, sum(sel)), col=colvec, pch=20, cex=.2)
				segments(dat$pos[ sel ] - .5, rep(e, sum(sel)), dat$pos[ sel ] + .5, rep(e, sum(sel)), col=colvec, pch=20, cex=.2)
		}
		}

		par(xpd=NA)	
		legend(max(d$pos[allsel]) + sum(allsel) * 0.05, ampMax, legend=c("Filtered nt", "Ref nt", "Non-ref nt"), col=c("#00000040", "#00000080", "red"), pch=15)
		legend(max(d$pos[allsel]) + sum(allsel) * 0.05, 0, legend=c("A", "C", "G", "T"), col=c("green", "blue", "black", "red"), pch=15)
		par(xpd=F)	
		par(mar=c(5, 4, 0, 7) + .1, xaxt="s")
		plot(NULL, type="n", xlab="Position", ylab=NA, xlim=range(d$pos[allsel]), ylim=c(0,4.5), yaxt="n")
		axis(side=2, at=1:4, labels=c("ins", "del", "snv", "ref"), las=2)
		abline(h=1:4, col="#00000040")
		sres <- nchar(as.character(resTmp$REF)) == nchar(as.character(resTmp$ALT))
		dres <- nchar(resTmp$REF) > nchar(resTmp$ALT)
		ires <- nchar(resTmp$REF) < nchar(resTmp$ALT)
		
		resTmp$marks <- rep(3.4, nrow(resTmp))
		
		resTmp$cols <- rep("#00000040", nrow(resTmp))
		resTmp$cols[grepl("^a", resTmp$ALT, ignore.case=T) & sres] <- "green"
		resTmp$cols[grepl("^c", resTmp$ALT, ignore.case=T) & sres] <- "blue"
		resTmp$cols[grepl("^g", resTmp$ALT, ignore.case=T) & sres] <- "black"
		resTmp$cols[grepl("^t", resTmp$ALT, ignore.case=T) & sres] <- "red"

		#resTmp$cols[dres] <- 2
		#resTmp$cols[ires] <- 2
		resTmp$marks[dres] <- 2.4
		resTmp$marks[ires] <- 1.4
		resTmp$markP1 <- resTmp$POS -.5
		resTmp$markP2 <- resTmp$POS +.5
		resTmp$markP1[dres] <- resTmp$POS[dres] +.5
		resTmp$markP2[dres] <- resTmp$POS[dres] +.5 + nchar(resTmp$REF[dres]) - 1
		resTmp$markP1[ires] <- resTmp$POS[ires] +.5
		resTmp$markP2[ires] <- resTmp$POS[ires] +.5 + 1
		ntCol <- rep(1, sum(allsel))
		ntCol[grepl("a", d$ntRef[allsel], ignore.case=T)] <- "green"
		ntCol[grepl("c", d$ntRef[allsel], ignore.case=T)] <- "blue"
		ntCol[grepl("g", d$ntRef[allsel], ignore.case=T)] <- "black"
		ntCol[grepl("t", d$ntRef[allsel], ignore.case=T)] <- "red"
	
		rect(d$pos[allsel] -.5, 4.4, d$pos[allsel] +.5, 3.6, col=ntCol, border=0)
		rect(resTmp$markP1, resTmp$marks, resTmp$markP2, resTmp$marks - .8, col=resTmp$cols, border=0)	
		
		if(length(clinPos != 0)) {
			rect(clinPos - .5, 3.5, clinPos + .5, .5, border="black", lwd=2)
		}


		dev.off()

#		png(paste(base, i, "vaf", "png", sep="."), width=960, height=480)
#		plot(d[sel & csel,2], d$vaf[sel & csel], type="l", ylim=c(0,100), main=paste(i, creg, cutoffs[3], sep=" "), xlab="position", ylab="%", col=2)
		#lines(d[sel & csel,2], d$pctVAR[sel & csel], type="l", col=2)
#		abline(h=100)
#		abline(h=c(1,5,10,50,90,95,99), lty=3, col=1)
#		abline(h=cutoffs[3], lty=4, col=2)
#		abline(v=c(limP[1], limP[1] + 150, limP[2] - 150, limP[2]),lty=4)
#		abline(v=c(qc[qcsel,8:9]), lty=5, col=2)
#		if (length(labs != 0)) {
#			text(resTmp$POS[resTmp$vaf > (cutoffs[3]/100)], 50, labels=resTmp$Amino_Acid_change[resTmp$vaf > (cutoffs[3]/100)])
#			abline(v=resTmp$POS[resTmp$vaf > (cutoffs[3]/100)], col="grey", lty=4)
#		}
#		dev.off()
	}
}

png(paste(base, "biasheat", "png", sep="."), width=960, height=960)
write.table(biasTable, paste(base, "bias", "tsv", sep="."), sep="\t")
biasTable[is.na(biasTable)] <- 0
heatmap(data.matrix(biasTable), Colv=NA, scale="r")
dev.off()
png(paste(base, "bias", "png", sep="."), width=960, height=480)
boxplot(data.matrix(biasTable), pch=20, cex=.5)
dev.off()
warnings()
#hist(d$pctVAR[d$pctVAR < 20], breaks=40, ylim=c(0,1000))
#hist(d$pctREF[d$pctREF > 80], breaks=40, ylim=c(0,1000))

#print(summary(d$pctVAR))


