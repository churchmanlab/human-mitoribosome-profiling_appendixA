

### Written 3/2024 by M. Couvillion


library(data.table)
library(stringr)
library(purrr)
library(rlist)
library('scales')
# library(gtsummary)

args <- commandArgs(trailingOnly = TRUE)

# path='/Users/Mary/Desktop/Data/hMitoRP/2024_BriBibel/hMitoRP/'
# Lib='De'
# sizerange='33to39'
# version='rpm'
# window=21
# use='fold'


path=paste0(getwd(),'/')
Lib = args[1]
sizerange = args[2] # 33to40 (std) 30to34 (sel)
version = args[3]  # rpm rps
window = as.numeric(args[4]) # 21
use = args[5]
suffixP = paste0('Asite_',sizerange,'_CodonCoverageEach_P.',version) # Asite_30to34_CodonCoverageEach_P.rps
suffixM = paste0('Asite_',sizerange,'_CodonCoverageEach_M.',version) # Asite_30to34_CodonCoverageEach_P.rps



sampDTP <- fread(paste0(path,'bedGraphs/',Lib,'_',suffixP, '.bedGraph'), skip=1)
sampDTM <- fread(paste0(path,'bedGraphs/',Lib,'_',suffixM, '.bedGraph'), skip=1)

# # Get values over sliding/rolling window
# plus
sampDT = sampDTP
slide.Mean = frollmean(sampDT$V4, window, algo='exact', align='center', na.rm=TRUE)
f = function(x) sd(x)
slide.SD = frollapply(sampDT$V4, window, f, align='center')
if (use == 'fold') {pausescore = sampDT$V4/slide.Mean}
if (use == 'SD') {pausescore = (sampDT$V4-slide.Mean)/slide.SD}

# # Make new DT with info above
sampDTP_2=sampDT[,.(V1=V1, V2=V2, V3=V3,V4=V4, slide.Mean=slide.Mean, slide.SD=slide.SD, pausescore=pausescore)]
sampDTP_2[is.na(sampDTP_2)] <- 0

# minus
# Make bookend data tables to make sure all of ND6 is included (windows have to be full-sized)
DTa <- data.table(V1='chrM', V2=seq(14700, 14673, by=-3), V3=seq(14703, 14676, by=-3), V4=0.00)
DTb <- data.table(V1='chrM', V2=seq(14145,14118, by=-3), V3=seq(14148,14121, by=-3), V4=0.00)
# Combine bookends with sample table
sampDT = rbind(DTa, sampDTM, DTb)
slide.Mean = frollmean(sampDT$V4, window, algo='exact', align='center', na.rm=TRUE)
f = function(x) sd(x)
slide.SD = frollapply(sampDT$V4, window, f, align='center')
if (use == 'fold') {pausescore = sampDT$V4/slide.Mean}
if (use == 'SD') {pausescore = (sampDT$V4-slide.Mean)/slide.SD}

# # Make new DT with info above
sampDTM_2=sampDT[,.(V1=V1, V2=V2, V3=V3,V4=V4, slide.Mean=slide.Mean, slide.SD=slide.SD, pausescore=pausescore)]
sampDTM_2[is.na(sampDTM_2)] <- 0


# Write bedGraphs


# plus
outpath = paste0(path, 'bedGraphs/', Lib,'_', suffixP,'_pauseScore',use,'_win',window,'.bedGraph')

cat("track type=bedGraph\n", file=outpath)
cols=c('V1', 'V2', 'V3', 'pausescore')
write.table(sampDTP_2[, ..cols], file=outpath, sep=("\t"), quote=FALSE, row.names=FALSE, append=TRUE,col.names=FALSE)

# minus
outpath = paste0(path, 'bedGraphs/', Lib,'_', suffixM,'_pauseScore',use,'_win',window,'.bedGraph')

cat("track type=bedGraph\n", file=outpath)
cols=c('V1', 'V2', 'V3', 'pausescore')
write.table(sampDTM_2[, ..cols], file=outpath, sep=("\t"), quote=FALSE, row.names=FALSE, append=TRUE,col.names=FALSE)

# source('~/Desktop/Data/hMitoRP/PanAnalysis/ScriptsForReanalysis/SlidingWindowPauseCaller.R')


