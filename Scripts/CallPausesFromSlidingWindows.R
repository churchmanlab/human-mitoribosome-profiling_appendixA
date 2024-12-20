library(data.table)


args <- commandArgs(trailingOnly = TRUE)


path=paste0(getwd(),'/')
Experiment= args[1] # CHLtest
soi = unlist(strsplit(args[2],","))
soiName = args[3]
sizeRange = args[4]
threshold= as.numeric(args[5]) # 2 for fold
use = args[6]
length = as.numeric(args[7])

# Experiment="CHLtest" # CHLtest
# soiName="Drug" # Drug
# threshold=2 # 2 for fold
# use="fold"

repnum=length(soi)

# Get data from bedgraphs
for (i in c(1:repnum)) {
assign(paste0('pDT',soi[i]), data.table(read.table(paste0(path,'bedGraphs/', soi[i], '_Asite_', sizeRange,'_CodonCoverageEach_P.rpm_pauseScore',use,'_win21.bedGraph'),skip=1)))
assign(paste0('mDT',soi[i]), data.table(read.table(paste0(path,'bedGraphs/', soi[i], '_Asite_', sizeRange,'_CodonCoverageEach_M.rpm_pauseScore',use,'_win21.bedGraph'),skip=1)))
}

# Combine all replicates
if (repnum == 1) {
pDT=cbind(get(paste0('pDT', soi[1])))
mDT=cbind(get(paste0('mDT', soi[1])))
colnames(pDT) <- c('chr', 'start', 'end', 's1')
colnames(mDT) <- c('chr', 'start', 'end', 's1')
}
if (repnum == 2) {
pDT=cbind(get(paste0('pDT', soi[1])), get(paste0('pDT', soi[2]))$V4)
mDT=cbind(get(paste0('mDT', soi[1])), get(paste0('mDT', soi[2]))$V4)
colnames(pDT) <- c('chr', 'start', 'end', 's1', 's2')
colnames(mDT) <- c('chr', 'start', 'end', 's1', 's2')
}
if (repnum == 3) {
pDT=cbind(get(paste0('pDT', soi[1])), get(paste0('pDT', soi[2]))$V4, get(paste0('pDT', soi[3]))$V4)
mDT=cbind(get(paste0('mDT', soi[1])), get(paste0('mDT', soi[2]))$V4, get(paste0('mDT', soi[3]))$V4)
colnames(pDT) <- c('chr', 'start', 'end', 's1', 's2', 's3')
colnames(mDT) <- c('chr', 'start', 'end', 's1', 's2', 's3')
}
if (repnum == 4) {
pDT=cbind(get(paste0('pDT', soi[1])), get(paste0('pDT', soi[2]))$V4, get(paste0('pDT', soi[3]))$V4, get(paste0('pDT', soi[4]))$V4)
mDT=cbind(get(paste0('mDT', soi[1])), get(paste0('mDT', soi[2]))$V4, get(paste0('mDT', soi[3]))$V4, get(paste0('mDT', soi[4]))$V4)
colnames(pDT) <- c('chr', 'start', 'end', 's1', 's2', 's3', 's4')
colnames(mDT) <- c('chr', 'start', 'end', 's1', 's2', 's3', 's4')
}
if (repnum == 5) {
pDT=cbind(get(paste0('pDT', soi[1])), get(paste0('pDT', soi[2]))$V4, get(paste0('pDT', soi[3]))$V4, get(paste0('pDT', soi[4]))$V4, get(paste0('pDT', soi[5]))$V4)
mDT=cbind(get(paste0('mDT', soi[1])), get(paste0('mDT', soi[2]))$V4, get(paste0('mDT', soi[3]))$V4, get(paste0('mDT', soi[4]))$V4, get(paste0('mDT', soi[5]))$V4)
colnames(pDT) <- c('chr', 'start', 'end', 's1', 's2', 's3', 's4', 's5')
colnames(mDT) <- c('chr', 'start', 'end', 's1', 's2', 's3', 's4', 's5')
}
if (repnum == 6) {
pDT=cbind(get(paste0('pDT', soi[1])), get(paste0('pDT', soi[2]))$V4, get(paste0('pDT', soi[3]))$V4, get(paste0('pDT', soi[4]))$V4, get(paste0('pDT', soi[5]))$V4, get(paste0('pDT', soi[6]))$V4)
mDT=cbind(get(paste0('mDT', soi[1])), get(paste0('mDT', soi[2]))$V4, get(paste0('mDT', soi[3]))$V4, get(paste0('mDT', soi[4]))$V4, get(paste0('mDT', soi[5]))$V4, get(paste0('mDT', soi[6]))$V4)
colnames(pDT) <- c('chr', 'start', 'end', 's1', 's2', 's3', 's4', 's5', 's6')
colnames(mDT) <- c('chr', 'start', 'end', 's1', 's2', 's3', 's4', 's5', 's6')
}




# Make column showing consensus (all above threshold)
up=threshold

if (repnum == 1) {
pDT[, pause := 0][s1 > up, pause := 10]
mDT[, pause := 0][s1 > up, pause := 10]
coi=c('s1')
}
if (repnum == 2) {
pDT[, pause := 0][s1 > up & s2 > up, pause := 10]
mDT[, pause := 0][s1 > up & s2 > up, pause := 10]
coi=c('s1','s2')
}
if (repnum == 3) {
pDT[, pause := 0][s1 > up & s2 > up & s3 > up, pause := 10]
mDT[, pause := 0][s1 > up & s2 > up & s3 > up, pause := 10]
coi=c('s1','s2', 's3')
}
if (repnum == 4) {
pDT[, pause := 0][s1 > up & s2 > up & s3 > up & s4 > up, pause := 10]
mDT[, pause := 0][s1 > up & s2 > up & s3 > up & s4 > up, pause := 10]
coi=c('s1','s2', 's3','s4')
}
if (repnum == 5) {
pDT[, pause := 0][s1 > up & s2 > up & s3 > up & s4 > up & s5 > up, pause := 10]
mDT[, pause := 0][s1 > up & s2 > up & s3 > up & s4 > up & s5 > up, pause := 10]
coi=c('s1','s2', 's3','s4','s5')
}
if (repnum == 6) {
pDT[, pause := 0][s1 > up & s2 > up & s3 > up & s4 > up & s5 > up & s6 > up, pause := 10]
mDT[, pause := 0][s1 > up & s2 > up & s3 > up & s4 > up & s5 > up & s6 > up, pause := 10]
coi=c('s1','s2', 's3','s4','s5', 's6')
}


# Add info for writing to file later
pDT$strand <- '+'
mDT$strand <- '-'

DT <- rbind(pDT, mDT)
DT <- DT[order(start)]

DT$name <- seq(1,nrow(DT))
DT$score <- rowMeans(DT[, ..coi])



# Get pauses
pause_pDT <- pDT[pause == 10]
pause_pDT$seqstart <- pause_pDT$start-length
pause_pDT$seqend <- pause_pDT$end+length
pause_pDT$strand <- '+'

pause_mDT <- mDT[pause == 10]
pause_mDT$seqstart <- pause_mDT$start-length
pause_mDT$seqend <- pause_mDT$end+length
pause_mDT$strand <- '-'

pause_DT <- rbind(pause_pDT, pause_mDT)
pause_DT <- pause_DT[order(start)]

pause_DT$name <- seq(1,nrow(pause_DT))
pause_DT$score <- '0'

print(paste0(soiName, ' ', use,' threshold ',threshold)) 
print(nrow(pause_DT))


bgoutfile=paste0(path,'bedGraphs/', soiName,'_slidingWindowPauses',use,'_Thresh',up,'_byCodon.bedGraph')
write('track type=bedGraph', bgoutfile)

# Write bedGraph for IGV
write.table(DT[, c('chr', 'start','end','pause')],file=bgoutfile, row.names=FALSE, quote=FALSE, sep='\t', col.names=FALSE, append=TRUE)


# Write bed for sample of interest for pulling out fasta sequence
boutfile=paste0(path, 'Pauses/',soiName,'_slidingWindowPauses',use,'_Thresh',up,'_byCodon.bed')
write.table(pause_DT[, c('chr', 'seqstart','seqend','name', 'score','strand')],file=boutfile, row.names=FALSE, quote=FALSE, sep='\t', col.names=FALSE, append=FALSE)


# source('/Users/Mary/Desktop/Data/hMitoRP/PanAnalysis/ScriptsForReanalysis/CallPausesFromSlidingWindows.R')