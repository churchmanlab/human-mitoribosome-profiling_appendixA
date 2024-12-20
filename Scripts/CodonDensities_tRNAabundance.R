library(data.table)
library('scales')
library(Rfast) # to use rowMins

# Experiment = 'CHLtest' # BBHEKsus BBHEKad GalShift MyoDiff1 OXPHOSinh HEK HEK1 MyoDiff GyanMitoRP1std
# site='E'
# soiName="Drug"
# n=4 # cutoff for codon frequency in genome


args <- commandArgs(trailingOnly = TRUE)

path=paste0(getwd(),'/')


Experiment = args[1] # 'CHLtest'
site = args[2] # 'EP'
soiName = args[3] # 'Drug'
n = as.numeric(args[4]) # 4


# Codon densities
if (Experiment == 'CHLtest') {
sizeRange = '32to38'
ctl = c("Hela_nodrug")
ctlName = 'Nodrug'
if (soiName == 'Drug') {
soi = c("Hela_drug")}
cdDT <- data.table(read.table(paste0(path, Experiment, '_', sizeRange,'_',site,'siteCodonDensities.txt'),header=TRUE))
ylimits = c(0,6)
}
if (Experiment == 'BBHEKsus') {
sizeRange = '33to39'
ctl = c("Da", "Db", "De", "Df")
ctlName = 'DMSO'
if (soiName == 'CHL') {
soi = c("Ca", "Cb", "Cc", "Ce", "Cf")}
if (soiName == 'LZD') {
soi = c("La", "Lb", "Lc", "Ld", "Le", "Lf")}
if (soiName == 'RET') {
soi = c("Ra", "Rb")}
cdDT <- data.table(read.table(paste0(path, Experiment, '_', sizeRange,'_',site,'siteCodonDensities.txt'),header=TRUE))
ylimits = c(0,6)
}

# For tRNA abundances and more codon info
taDT <- data.table(read.table(paste('./SeqFiles/hMitoGeneticCode_BattersbyMyoblast_tRNAabundance.txt'), header=TRUE, , sep = ','))

# Merge
DT = data.table(merge(cdDT,taDT, by='codon'))
# Keep only codons with n or more instances in transcriptome
DT = DT[DT$CodonFreq >= n]

setorder(DT, cols = 'AminoAcid')

SampNum = (ncol(DT) - 7) / 2

dens = paste0(soi, '_aveDensity')
sd = paste0(soi, '_SD')

ctldens = paste0(ctl, '_aveDensity')
ctlsd = paste0(ctl, '_SD')

means = rowMeans(DT[, ..dens])
ranges = rowMaxs(as.matrix(DT[, ..dens]), value=TRUE) - rowMins(as.matrix(DT[, ..dens]), value=TRUE)

ctlmeans = rowMeans(DT[, ..ctldens])
ctlranges = rowMaxs(as.matrix(DT[, ..ctldens]), value=TRUE) - rowMins(as.matrix(DT[, ..ctldens]), value=TRUE)

# plotDT = DT[, ..dens]
# sdDT = DT[, ..sd]


colors = c('red', 'darkorange','darkgoldenrod1', 'forestgreen','dodgerblue', 'blue','darkblue','purple','violet','grey50','black','brown','red', 'darkorange','darkgoldenrod1', 'forestgreen','dodgerblue', 'darkblue','purple','violet') #,'grey50','black')

A=nrow(DT[grep('^A', oneLetter)])
C=nrow(DT[grep('^C', oneLetter)])
D=nrow(DT[grep('^D', oneLetter)])
E=nrow(DT[grep('^E', oneLetter)])
F=nrow(DT[grep('^F', oneLetter)])
G=nrow(DT[grep('^G', oneLetter)])
H=nrow(DT[grep('^H', oneLetter)])
I=nrow(DT[grep('^I', oneLetter)])
K=nrow(DT[grep('^K', oneLetter)])
L=nrow(DT[grep('^L', oneLetter)])
M=nrow(DT[grep('^M', oneLetter)])
N=nrow(DT[grep('^N', oneLetter)])
P=nrow(DT[grep('^P', oneLetter)])
Q=nrow(DT[grep('^Q', oneLetter)])
R=nrow(DT[grep('^R', oneLetter)])
S=nrow(DT[grep('^S', oneLetter)])
T=nrow(DT[grep('^T', oneLetter)])
V=nrow(DT[grep('^V', oneLetter)])
W=nrow(DT[grep('^W', oneLetter)])
Y=nrow(DT[grep('^Y', oneLetter)])
repTimes = c(A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V) # * length(soi)

names = paste0(DT$codon,' ',DT$AminoAcid)

################ Barplot of codon densities ################
pdf(paste0(path,'plots/', Experiment,'_', soiName, '_', sizeRange,'_',site,'siteCodonDens_barplot.pdf'),
     width=6.8, # width=4.5, width=2, width=2.5
     height=3, #4
     pointsize=11)

xx = barplot(ctlmeans, beside=TRUE, names = names, col = rep(colors,repTimes), border=NA, ylab = 'Mean codon occupancy', cex.axis=.7, cex.names = .65, ylim=ylimits, las=2, main=paste0(ctlName,'\n',site,' site')) # t(plotDT)

if (length(ctl) > 1) {
# range bars
arrows(xx, ctlmeans-(ctlranges/2), xx, ctlmeans+(ctlranges/2), length=0.03, angle=90, code=3, lwd=.2)
}

xx = barplot(means, beside=TRUE, names = names, col = rep(colors,repTimes), border=NA, ylab = 'Mean codon occupancy', cex.axis=.7, cex.names = .65, ylim=ylimits, las=2, main = paste0(soiName,'\n',site,' site')) # t(plotDT)

if (length(soi) > 1) {
# range bars
arrows(xx, means-(ranges/2), xx, means+(ranges/2), length=0.03, angle=90, code=3, lwd=.2)
}

# standard deviation
# arrows(xx, t(plotDT)-t(sdDT), xx, t(plotDT)+t(sdDT), length=0.01, angle=90, code=3, lwd=.2)

abline(h=1, lwd = .2, lty = 2, col = alpha('grey60', .5))

# Plot SD

# x=c(1.5:((nrow(DT)*(SampNum+1))+.5))
# # Make 0 column at end of data table to make error bars line up
# fill_plotDT <- copy(plotDT)
# fill_plotDT[, fill := 0]
# y=as.vector(t(fill_plotDT))
# # Get sd
# sdIdx = seq(3,2*SampNum+2,2)
# sdDT <- DT[,sdIdx, with=FALSE]
# # Make 0 column at end of data table to make error bars line up
# fill_sdDT <- copy(sdDT)
# fill_sdDT[, fill := 0]
# sd=as.vector(t(fill_sdDT))
# 
# arrows(x, y-sd, x, y+sd, length=0.005, angle=90, code=3, lwd = 0.1)



################ Barplot of fold-change ################

xx = barplot(log(means/ctlmeans,2), beside=TRUE, names = names, col = rep(colors,repTimes), main = paste0(soiName,'\n',site,' site'), border=NA, ylab = paste0('Mean codon occupancy ',soiName,'/',ctlName, ' (log2)'), cex.axis=.7, cex.names = .65, las=2, ylim=c(-1,2)) # t(plotDT)

abline(h=1, lwd = .2, lty = 2, col = alpha('grey60', .5))


dev.off()





################ Individual codon densities ################
# libName = soi[1]
# indDT <- data.table(read.table(paste0(path, libName, '_', sizeRange,'_indCodonDensities.txt'),header=FALSE, sep='['))
# 
# ## Asp
# codon = 'GAT'
# occ = as.character(indDT[V1==paste0(codon,',')]$V2)
# GAToccs = as.double(strsplit(strsplit(occ,']')[[1]],', ')[[1]]) # Split by ',' and ']'
# codon = 'GAC'
# occ = as.character(indDT[V1==paste0(codon,',')]$V2)
# GACoccs = as.double(strsplit(strsplit(occ,']')[[1]],', ')[[1]]) # Split by ',' and ']'
# Asp_occs = c(GAToccs, GACoccs)
# 
# # Leu
# codon = 'TTG'
# occ = as.character(indDT[V1==paste0(codon,',')]$V2)
# TTGoccs = as.double(strsplit(strsplit(occ,']')[[1]],', ')[[1]]) # Split by ',' and ']'
# codon = 'TTA'
# occ = as.character(indDT[V1==paste0(codon,',')]$V2)
# TTAoccs = as.double(strsplit(strsplit(occ,']')[[1]],', ')[[1]]) # Split by ',' and ']'
# Leu_occs = c(TTGoccs, TTAoccs)
# 
# # Phe
# codon = 'TTT'
# occ = as.character(indDT[V1==paste0(codon,',')]$V2)
# TTToccs = as.double(strsplit(strsplit(occ,']')[[1]],', ')[[1]]) # Split by ',' and ']'
# codon = 'TTC'
# occ = as.character(indDT[V1==paste0(codon,',')]$V2)
# TTCoccs = as.double(strsplit(strsplit(occ,']')[[1]],', ')[[1]]) # Split by ',' and ']'
# Phe_occs = c(TTToccs, TTCoccs)
# 
# # Arg
# codon = 'CGT'
# occ = as.character(indDT[V1==paste0(codon,',')]$V2)
# CGToccs = as.double(strsplit(strsplit(occ,']')[[1]],', ')[[1]]) # Split by ',' and ']'
# Arg_occs = CGToccs
# 
# # His
# codon = 'CAC'
# occ = as.character(indDT[V1==paste0(codon,',')]$V2)
# CACoccs = as.double(strsplit(strsplit(occ,']')[[1]],', ')[[1]]) # Split by ',' and ']'
# His_occs = CACoccs
# 
# 
# pdf(paste0(path,'plots/', Experiment, '_', sizeRange,'_',site,'siteindCodonDens_boxplot.pdf'),
#      width=4, # width=4.5, width=2, width=2.5
#      height=4, #4
#      pointsize=11)
#      
# boxplot(Asp_occs, Leu_occs, Phe_occs, Arg_occs, His_occs,names = c('Asp', 'Leu', 'Phe', 'Arg', 'His'), outline = FALSE, ylab = 'Codon Occupancy', ylim = c(0,16)) # ylim = c(0,9)
# points(jitter(rep(1,length(Asp_occs)),12), Asp_occs, pch = 16,cex = .7, col = alpha('forestgreen',.6))
# points(jitter(rep(2,length(Leu_occs)),5), Leu_occs, pch = 16,cex = .7, col = alpha('black', .6))
# points(jitter(rep(3,length(Phe_occs)),3), Phe_occs, pch = 16,cex = .7, col = alpha('darkorange', .6))
# points(jitter(rep(4,length(Arg_occs)),3), Arg_occs, pch = 16,cex = .7, col = alpha('darkorange', .6))
# points(jitter(rep(5,length(His_occs)),3), His_occs, pch = 16,cex = .7, col = alpha('violet', .6))
# 
# abline(h=1, lty = 2, lwd = 0.2)
# 
# dev.off()
# 
# Barplot of codon occupancies vs tRNA abundance
pdf(paste0(path,'plots/', Experiment,'_', soiName, '_', sizeRange,'_',site,'siteCodonDens_vs_tRNAabundance.pdf'),
     width=4.5, # width=4.5, width=2, width=2.5
     height=3.5, #4
     pointsize=11)

tRNAabund = DT$meanAbundance
codonOcc = means

pearson_r = cor.test(~tRNAabund[tRNAabund>0]+codonOcc[tRNAabund>0], method = c('pearson')) 
mylabel = bquote(italic(r) == .(format(pearson_r$estimate, digits = 2)))

plot(tRNAabund, codonOcc, main = soi, cex.main = 0.5, cex.axis = 1, col = rep(colors,repTimes), pch = 16, cex = 1, xlab = 'tRNA abundance (tpm)', ylab = 'Mean codon occupancy', ylim=ylimits) 


# tRNA abundance SD (3 replicates)
arrows(tRNAabund-DT$SD, codonOcc, tRNAabund+DT$SD, codonOcc, length=0.03, angle=90, code=3, lwd=.3)
# codon occupancy ranges
arrows(tRNAabund, codonOcc-(ranges/2), tRNAabund, codonOcc+(ranges/2), length=0.03, angle=90, code=3, lwd=.3)

# plot points again to make them on top of error bars
points(tRNAabund, codonOcc,col = rep(colors,repTimes), pch = 16, cex = 1.5)

text(x = .002, y = max(ylimits)-.5, labels = mylabel, cex = 1)

dev.off()

# source('/Users/Mary/Desktop/Data/hMitoRP/PanAnalysis/ScriptsForReanalysis/CodonDensities_tRNAabundance.R')

