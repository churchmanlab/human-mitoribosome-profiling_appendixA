library(data.table)
library('scales')
library(Rfast) # to use rowMins

# Experiment = 'CHLtest'
# site='EP'
# soiName='Drug'

args <- commandArgs(trailingOnly = TRUE)

path=paste0(getwd(),'/')


Experiment = args[1] # 'CHLtest'
site = args[2] # 'EP'
soiName = args[3] # 'Drug'
n = as.numeric(args[4]) # have to be present at least this many times in transcriptome

# Codon densities
if (Experiment == 'CHLtest') {
sizeRange = '32to38'
con = c("Hela_nodrug")
conName = 'Nodrug'
if (soiName == 'Drug') {
soi = c("Hela_drug")}
cdDT <- data.table(read.table(paste0(path, Experiment, '_', sizeRange,'_',site,'siteCodonDensities.txt'),header=TRUE))
ylimits = c(0,10)
}
if (Experiment == 'BBHEKsus') {
sizeRange = '33to39'
con = c("Da", "Db", "De", "Df")
conName = 'DMSO'
if (soiName == 'CHL') {
soi = c("Ca", "Cb", "Cc", "Ce", "Cf")}
if (soiName == 'LZD') {
soi = c("La", "Lb", "Lc", "Ld", "Le", "Lf")}
if (soiName == 'RET') {
soi = c("Ra", "Rb")}
cdDT <- data.table(read.table(paste0(path, Experiment, '_', sizeRange,'_',site,'siteCodonDensities.txt'),header=TRUE))
ylimits = c(0,10)
}



# For tRNA abundances and more codon info
taDT <- data.table(read.table(paste('./SeqFiles/hMitoGeneticCode_BattersbyMyoblast_tRNAabundance.txt'), header=TRUE, , sep = ','))
codeDT <- taDT[, c('codon', 'oneLetter')]

# Only keep those where codon frequency > n
DT=cdDT[cdDT$CodonFreq >= n]

DT[, codon1 := substring(codon,1,3)]
DT[, codon2 := substring(codon,4,6)]

# Look up first codon one-letter code and add to DT
setkey(DT, codon1)
setkey(codeDT, codon)
DT = DT[codeDT]
# Look up second codon one-letter code and add to DT
setkey(DT, codon2)
setkey(codeDT, codon)
DT = DT[codeDT]
# Remove NAs
DT <- na.omit(DT)
# Make column with dicodon one letter code
DT[, dicodon := paste0(oneLetter, i.oneLetter)]
# Remove NAs
DT <- na.omit(DT)

# Alphabetical order
setorder(DT, cols = 'dicodon')

dens = paste0(soi, '_aveDensity')
sd = paste0(soi, '_SD')
condens = paste0(con, '_aveDensity')
consd = paste0(con, '_SD')
# Keep only columns of interest
DT <- DT[,c(condens, dens, 'dicodon', 'CodonFreq'), with=FALSE]


# Sum Codon frequency, take mean of codon density by dicodon
freq <- DT[,.(CodonFreq.Sum=sum(CodonFreq)),by=dicodon]
DTagg <- DT[, lapply(.SD, mean), by= dicodon]
DTagg[, CodonFreq := freq$CodonFreq.Sum]
# Keep only dicodons that occur more than n times
DTagg <- DTagg[CodonFreq > 8]

dicodon=DTagg$dicodon

means = rowMeans(DTagg[, ..dens])
conmeans = rowMeans(DTagg[, ..condens])

ranges = rowMaxs(as.matrix(DTagg[, ..dens]), value=TRUE) - rowMins(as.matrix(DTagg[, ..dens]), value=TRUE)
conranges = rowMaxs(as.matrix(DTagg[, ..condens]), value=TRUE) - rowMins(as.matrix(DTagg[, ..condens]), value=TRUE)


fold = means/conmeans

pdf(paste0(path,'plots/', Experiment,'_', soiName, '_', sizeRange,'_',site,'siteCodonDens_foldchange_barplot.pdf'),
     width=15, # width=4.5, width=2, width=2.5
     height=4, #4
     pointsize=11)

colors = c('red', 'darkorange','darkgoldenrod1', 'forestgreen','dodgerblue', 'blue','darkblue','purple','violet','grey50','black','brown','red', 'darkorange','darkgoldenrod1', 'forestgreen','dodgerblue', 'darkblue','purple','violet')


A=nrow(DTagg[grep('^A', dicodon)])
C=nrow(DTagg[grep('^C', dicodon)])
D=nrow(DTagg[grep('^D', dicodon)])
E=nrow(DTagg[grep('^E', dicodon)])
F=nrow(DTagg[grep('^F', dicodon)])
G=nrow(DTagg[grep('^G', dicodon)])
H=nrow(DTagg[grep('^H', dicodon)])
I=nrow(DTagg[grep('^I', dicodon)])
K=nrow(DTagg[grep('^K', dicodon)])
L=nrow(DTagg[grep('^L', dicodon)])
M=nrow(DTagg[grep('^M', dicodon)])
N=nrow(DTagg[grep('^N', dicodon)])
P=nrow(DTagg[grep('^P', dicodon)])
Q=nrow(DTagg[grep('^Q', dicodon)])
R=nrow(DTagg[grep('^R', dicodon)])
S=nrow(DTagg[grep('^S', dicodon)])
T=nrow(DTagg[grep('^T', dicodon)])
V=nrow(DTagg[grep('^V', dicodon)])
W=nrow(DTagg[grep('^W', dicodon)])
Y=nrow(DTagg[grep('^Y', dicodon)])
repTimes = c(A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y) # * length(soi)

xx = barplot(conmeans, beside=TRUE, names = dicodon, col = rep(colors,repTimes), border=NA, ylab = 'Mean dicodon occupancy', cex.axis=.7, cex.names = .65, ylim=ylimits, las=2, main=paste0(conName,'\n',site,' site')) # t(plotDT)

if (length(con) > 1) {
# range bars
arrows(xx, conmeans-(conranges/2), xx, conmeans+(conranges/2), length=0.03, angle=90, code=3, lwd=1)
}



xx = barplot(means, beside=TRUE, names = dicodon, col = rep(colors,repTimes), border=NA, ylab = 'Mean dicodon occupancy', cex.axis=.7, cex.names = .65, ylim=ylimits, las=2, main=paste0(soiName,'\n',site,' site')) # t(plotDT)

if (length(soi) > 1) {
# range bars
arrows(xx, means-(ranges/2), xx, means+(ranges/2), length=0.03, angle=90, code=3, lwd=1)
}

xx = barplot(log(fold,2), beside=TRUE, names = dicodon, col = rep(colors,repTimes), border=NA, ylab = paste0('Mean codon occupancy ',soiName,'/',conName, ' (log2)'), cex.axis=.7, cex.names = .65, ylim=c(-2.1,2.1), las=2, main=paste0(soiName,'\n',site,' site')) # t(plotDT)




dev.off()



# source('/Users/Mary/Desktop/Data/hMitoRP/PanAnalysis/ScriptsForReanalysis/diCodonDensities_tRNAabundance_foldchange.R')

