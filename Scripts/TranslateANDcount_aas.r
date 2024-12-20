library(data.table)
library(Biostrings)
library(matrixStats)
library(scales)



args <- commandArgs(trailingOnly = TRUE)

type = args[1] #'slidingWindowPauses' # reads slidingWindowPauses


colors=c('royalblue','gold1','indianred2','seagreen3')



# old tables:
# e.g. P-5 to P+1 is c(-6,0)
# -5 -4 -3 -2 -1 0 1
#  6  5  4  3  2 1 0

# new tables:
# 15 amino acids in the table
# -6   -5  -4  -3  -2 -1  P  A  1  2  3  4  5  6 7
# -14 -13 -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0

Ppositions=c('P-6','P-5','P-4','P-3', 'P-2', 'P-1', 'P', 'P+1', 'P+2','P+3','P+4','P+5','P+6','P+7','P+8')
# Initialize vector for keeping track of info content at each position
Cinfos=c()
Linfos=c()


pos=c(14:0)

for (n in pos) {
	positions = c(-n,-n)
	
	P=Ppositions[(length(Ppositions)-n)]


	if (type=='reads') {
	soiNames=c('Da','Db','De', 'Df','Ca', 'Cb', 'Cc','Ce','Cf','La', 'Lb', 'Ld','Le','Lf', 'Ra', 'Rb', 'Dg','Dh','Ia','Ib')
	path='./bedGraphs/fastabeds/'
	}
	if (type=='pauses') {
	soiNames=c('CHL','LZD','IND')
	path='./Pauses/'
	}
	if (type=='slidingWindowPauses') {
	soiNames=c('DMSO','CHL','LZD', 'IND')
	path='./Pauses/'	
	}



	freqss=list()
	successess=list() # To keep track of the numbers for the Two Proportions test
	totals=c()
	for (i in c(1:length(soiNames))) {

		soiName=soiNames[i]

		if (soiName == 'Dg' | soiName == 'Dh' | soiName == 'Ia' | soiName == 'Ib' | soiName == 'IND') {
		sizeRange = '29to36'} else {
		sizeRange = '33to39'}


		if (type=='reads') {
		DT <- fread(paste0(path, soiName,'_Asite_',sizeRange,'_CodonCoverageEach.rpm.fasta'))
		}
		if (type=='pauses') {
		DT <- fread(paste0(path, soiName,'_differencePauses_byCodon.tab.fasta'))
		}
		if (type=='slidingWindowPauses') {
		DT <- fread(paste0(path, soiName,'_slidingWindowPausesfold_Thresh2_byCodon.tab.fasta'))
		}

		# Translate
		SGC1 <- getGeneticCode("SGC1")
		seqs=DT$V7
		pep = as.data.table(translate(DNAStringSet(DT$V7), genetic.code=SGC1, no.init.codon=TRUE))$x
		DT[, pep := pep ]


		# Write table for counting amino acids
		if (type=='reads') {

		# Keep only positive values
		DT <- DT[DT$V5 > 0]
		write.table(DT[, c('pep', 'V5')],file=paste0(path, soiName,'_Asite_',sizeRange,'_CodonCoverageEach.rpm.pep.tab'), row.names=FALSE, quote=FALSE, sep='\t', col.names=FALSE, append=FALSE)
		}
		if (type=='pauses') {
		write.table(DT[, c('pep', 'V5')],file=paste0(path, soiName,'_differencePauses_byCodon.pep.tab'), row.names=FALSE, quote=FALSE, sep='\t', col.names=FALSE, append=FALSE)
		}
		if (type=='slidingWindowPauses') {
		write.table(DT[, c('pep', 'V5')],file=paste0(path, soiName,'_slidingWindowPausesfold_Thresh2_byCodon.pep.tab'), row.names=FALSE, quote=FALSE, sep='\t', col.names=FALSE, append=FALSE)
		}

		# Count each amino acid in all ribosome positions listed

		aas = c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V')

		if (type=='reads') {
			total=sum(DT$V5)
			totals=c(totals, total)

			successes=c()
			freqs=c()
			for (aa in aas) {
				success=sum(DT[grep(aa,substr(DT$pep, nchar(DT$pep)+positions[1], nchar(DT$pep)+positions[2]))]$V5)
				successes=c(successes, success)
				freq=success/total*100
				freqs=c(freqs,freq)
			
			}
		}

		if (type=='pauses' | type=='slidingWindowPauses') {
			total=nrow(DT)
			totals=c(totals, total)


			successes=c()
			freqs=c()
			for (aa in aas) {
			success=nrow(DT[grep(aa,substr(DT$pep, nchar(DT$pep)+positions[1], nchar(DT$pep)+positions[2]))])
			freq=success/total*100
			successes=c(successes, success)
			freqs=c(freqs,freq)
			}
		}

		freqss[[length(freqss) + 1]] <- freqs
		successess[[length(successess) + 1]] <- successes
		} # Finish going through each sample/replicate



	pdf(paste0(path, 'plots/AAfrequency_',P,'_',type,'.pdf'), 
		 width=8, # width=4.5, width=2, width=2.5
		 height=3.5, #4
		 pointsize=11)
  
	if (type=='reads') {
		mat=matrix(unlist(freqss), nrow=length(aas))
		D=mat[,c(1:4)]
		Dmean=rowMeans(D)
		Dsd=rowSds(D)
		C=mat[,c(5:9)]
		Cmean=rowMeans(C)
		Csd=rowSds(C)
		L=mat[,c(10:14)]
		Lmean=rowMeans(L)
		Lsd=rowSds(L)
		
		# T-test for C vs D
		Cps=c()
		Lps=c()
		for (r in c(1:nrow(D))) {
			Cp=t.test(C[r,],D[r,])$p.value
			Lp=t.test(L[r,],D[r,])$p.value
			Cps = c(Cps, Cp)
			Lps = c(Lps, Lp)
		}
		Cpsigs=which(Cps<.005)
		Lpsigs=which(Lps<.005)
		
		
		
		toPlot=matrix(c(Dmean,Cmean, Lmean), nrow=length(aas))
		sds=matrix(c(Dsd,Csd, Lsd), nrow=length(aas))
		sucmat=matrix(unlist(successess), nrow=length(aas))


		title=paste0('Fraction of reads with each aa in ',P)
		legtext=c('DMSO','CHL','LZD')
		
	
		# Get sum of log2 fold changes
# 		Cinfo <- sum(abs(log(Cmean/Dmean, 2)))
# 		Linfo <- sum(abs(log(Lmean/Dmean, 2)))
# 		Cinfos <- c(Cinfos, Cinfo)
# 		Linfos <- c(Linfos, Linfo)
		
		Cinfo <- sum(pmax(Cmean, Dmean)-pmin(Cmean, Dmean))
		Linfo <- sum(pmax(Lmean, Dmean)-pmin(Lmean, Dmean))
		Cinfos <- c(Cinfos, Cinfo)
		Linfos <- c(Linfos, Linfo)
	
		# 2-sample test for equality of proportions with continuity correction or not (correct = FALSE)
		# https://www.r-bloggers.com/2022/05/two-sample-proportions-test-in-r-complete-guide/
		# https://www.statisticshowto.com/probability-and-statistics/hypothesis-testing/z-test/
		s1pvals=c()
		s2pvals=c()
		s3pvals=c()
		for (a in c(1:length(aas))) {
			res1 <- prop.test(x = c(sucmat[a,5:9], sucmat[a,1:4]), n = c(totals[5:9], totals[1:4]))
			res2 <- prop.test(x = c(sucmat[a,10:14], sucmat[a,1:4]), n = c(totals[5:9], totals[1:4]))
		
			s1pval=res1$p.value
			s2pval=res2$p.value
			s1pvals=c(s1pvals,s1pval)
			s2pvals=c(s2pvals,s2pval)
		
		}
		readsig1 = which(s1pvals<.0001)
		readsig2 = which(s2pvals<.0001)
	

	
	
	
	}


	if (type=='pauses' | type=='slidingWindowPauses') {
		mat=matrix(unlist(freqss), nrow=length(aas))
		D=mat[,c(1)]
		Dmean=D
		Dsd=0
		C=mat[,c(2)]
		Cmean=C
		Csd=0
		L=mat[,c(3)]
		Lmean=L
		Lsd=0
		if (length(soiNames) == 4) {
			I=mat[,c(4)]
			Imean=I
			Isd=0
		}


		toPlot=matrix(c(Dmean,Cmean, Lmean), nrow=length(aas))
		sds=matrix(c(Dsd,Csd, Lsd), nrow=length(aas)) # there is only one, so these will be 0
		if (length(soiNames) == 4) {
			toPlot=matrix(c(Dmean,Cmean, Lmean, Imean), nrow=length(aas))
			sds=matrix(c(Dsd,Csd, Lsd, Isd), nrow=length(aas))
		}
		if (type == 'slidingWindowPauses') {
			sucmat=matrix(unlist(successess), nrow=length(aas))
		}

		title=paste0('Fraction of pauses with each aa in ',P)
		legtext=soiNames
			if (type=='slidingWindowPauses') {
				if (length(soiNames) == 3) {
					legtext=c('DMSO','CHL','LZD')
				}
				if (length(soiNames) == 4) {
					legtext=c('DMSO','CHL','LZD','IND')
				}
			}
	}


	if (type == 'slidingWindowPauses') {
		# 2-sample test for equality of proportions with continuity correction
		# https://www.r-bloggers.com/2022/05/two-sample-proportions-test-in-r-complete-guide/
		# https://www.statisticshowto.com/probability-and-statistics/hypothesis-testing/z-test/
		s1pvals=c()
		s2pvals=c()
		s3pvals=c()
		for (i in c(1:length(aas))) {
			res1 <- prop.test(x = c(sucmat[i,1], sucmat[i,2]), n = c(totals[1], totals[2]))
			res2 <- prop.test(x = c(sucmat[i,1], sucmat[i,3]), n = c(totals[1], totals[3]))
			if (length(soiNames) == 4) {
				res3 <- prop.test(x = c(sucmat[i,1], sucmat[i,4]), n = c(totals[1], totals[4]))
			}
			s1pval=res1$p.value
			s2pval=res2$p.value
			s1pvals=c(s1pvals,s1pval)
			s2pvals=c(s2pvals,s2pval)
			if (length(soiNames) == 4) {
				s3pval=res3$p.value
				s3pvals=c(s3pvals,s3pval)
			}
		
		}
		pausesig1 = which(s1pvals<.01)
		pausesig2 = which(s2pvals<.01)
		if (length(soiNames) == 4) {
			pausesig3 = which(s3pvals<.01)
		}
	}


	xx = barplot(t(toPlot), beside=TRUE, names=aas, border=NA, col=colors[1:length(soiNames)], cex.axis=.8,main=title, legend.text=legtext,args.legend = list(x='topright', bty = 'n', fill = colors[1:length(soiNames)], cex=.9), ylim=c(0,36))

	# Bars for standard deviation across replicates for reads measure
	arrows(xx, as.vector(t(toPlot))-as.vector(t(sds))/2, xx, as.vector(t(toPlot))+as.vector(t(sds))/2, length=0.05, angle=90, code=3, lwd=.5)
	
	if (type == 'reads') {
		# Indicate significance for reads
		
		# Test of 2-proportions
# 		xx1 = xx[,readsig1, drop=FALSE][2,]
# 		yy1 = t(toPlot)[,readsig1, drop=FALSE][2,]
# 		xx2 = xx[,readsig2, drop=FALSE][3,]
# 		yy2 = t(toPlot)[,readsig2, drop=FALSE][3,]
# 		if (length(xx1) >0) {
# 			text(xx1, yy1+1, labels='*', col='black')
# 		}
# 		if (length(xx2) >0) {
# 			text(xx2, yy2+1, labels='*', col='black')
# 		}


		# T test
# 		xx1 = xx[2,Cpsigs]
# 		yy1 = t(toPlot)[2,Cpsigs]
# 		xx2 = xx[3,Lpsigs]
# 		yy2 = t(toPlot)[3,Lpsigs]
# 		if (length(xx1) > 0) {
# 			text(xx1, yy1+2, labels='*', col='black')
# 		}
# 		if (length(xx2) > 0) {
# 			text(xx2, yy2+2, labels='*', col='black')
# 		}
	}	

	if (type == 'slidingWindowPauses') {
		# Indicate significance for sliding window pauses
		xx1 = xx[,pausesig1, drop=FALSE][2,]
		yy1 = t(toPlot)[,pausesig1, drop=FALSE][2,]
		xx2 = xx[,pausesig2, drop=FALSE][3,]
		yy2 = t(toPlot)[,pausesig2, drop=FALSE][3,]
		if (length(xx1) >0) {
			text(xx1, yy1+1, labels='*', col='black')
		}
		if (length(xx2) >0) {
			text(xx2, yy2+1, labels='*', col='black')
		}
		if (length(soiNames) == 4) {
			xx3 = xx[,pausesig3, drop=FALSE][4,]
			yy3 = t(toPlot)[,pausesig3, drop=FALSE][4,]
			if (length(xx3) >0) {
				text(xx3, yy3+1, labels='*', col='black')
			}
		}
		
	}

	dev.off()

	toWrite=data.table(toPlot)
	# rownames(toWrite) <- soiNames
	# colnames(toWrite) <- aas
	rownames(toWrite) <- aas
	colnames(toWrite) <- legtext

	fwrite(toWrite,file=paste0(path, 'plots/AAfrequency_',P,'_',type,'.txt'), row.names=TRUE, quote=FALSE, sep='\t', col.names=TRUE)

}

if (type == 'reads') {
	pdf(paste0(path, 'plots/infocontent_reads.pdf'), 
		 width=4, # width=4.5, width=2, width=2.5
		 height=3.5, #4
		 pointsize=11)
		 
		 
	# Normalize to P+8 where there shouldn't be any effect (also try the minimum)
	CinfosNorm = abs(Cinfos-Cinfos[length(Cinfos)]) # min(Cinfos)
	LinfosNorm = abs(Linfos-Linfos[length(Linfos)]) #min(Linfos)
	# Now calculate fraction at each position
	CinfosFrac = CinfosNorm/sum(CinfosNorm)*100
	LinfosFrac = LinfosNorm/sum(LinfosNorm)*100
	
	toPlot=matrix(c(CinfosFrac, LinfosFrac), ncol=2)
	colvec=colors[c(2:3)]
	barplot(t(toPlot), beside=TRUE,col = colvec, border=NA, ylab = 'Fraction predictive information', names=Ppositions[c((length(Ppositions)-pos[1]):(length(Ppositions)-pos[length(pos)]))], las=2)
	legend("topleft", legend = c('CHL','LZD'), border='white',fill = colvec, bty = "n", cex = 0.7)


	dev.off()
}



# source('/Users/Mary/Desktop/Data/hMitoRP/2024_BriBibel/Scripts/TranslateANDcount_aas.r')
