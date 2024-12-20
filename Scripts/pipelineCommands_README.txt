




1. A site transform and count frame

### Input note ### 
### This script requires a bam file of aligned single-end reads named as ${lib}_Aligned.Mito_mRNA.noDups.bam. See the repository "human-mitoribosome-profiling" to determine appropriate sizes and offsets

### Output note ### 
### This script outputs many files including A-site bedgraphs, periodicity counts with ignoring either the first 3 or 6 codons, coverage counts on each codon for downstream analyses (*_CodonCoverageEach.txt), as well as a version to load on IGV (_CodonCoverageEach_P/M.bedGraph), and reads per kb

Exp="CHLtest"
Libs="Hela_nodrug Hela_drug"
sizerange='32to38'
sizes="32 33 34 35 36 37 38"
offsets="15 16 16 17 18 19 20"

for lib in $Libs
do
sbatch /n/groups/churchman/mc348/hMitoRP/PanAnalysis/Scripts/AsiteAndCountFrame.sh $lib $Exp $sizerange "$sizes" "$offsets"
done


# Combine outputs for periodicity counts

Experiment="CHLtest"
Libs=("Hela_nodrug" "Hela_drug")
sizeRange='32to38'

paste <(awk '{print $2"\t"$3}' ${Libs[0]}_Asite_${sizeRange}_FrameCount_ignore1st6.txt) <(awk '{print $3}' ${Libs[1]}_Asite_${sizeRange}_FrameCount_ignore1st6.txt) > ${Experiment}_Asite_${sizeRange}_FrameCount_ignore1st6.txt






2a. Normalize bedGraphs by rpm

Inputs: 
${LibName}_${suffix}_P.bedGraph
${LibName}_${suffix}_PAll.txt
${LibName}_${suffix}_M.bedGraph
${LibName}_${suffix}_MAll.txt
Outputs: 
${LibName}_${suffix}_P.rpm.bedGraph
${LibName}_${suffix}_PAll.rpm.txt
${LibName}_${suffix}_P.rps.bedGraph
${LibName}_${suffix}_PAll.rps.txt
${LibName}_${suffix}_M.rpm.bedGraph
${LibName}_${suffix}_MAll.rpm.txt
${LibName}_${suffix}_M.rps.bedGraph
${LibName}_${suffix}_MAll.rps.txt


Exp="CHLtest"
Libs=("Hela_nodrug" "Hela_drug")
sizeRange="32to38"
suffix=Mito_mRNA.noDups.Asite_${sizeRange}
spikes=(3000 3000) # If no spike ins here, just add this so it will run

for i in {0..1}
do
sbatch ./Scripts/NormalizeBedgraph.sh ${Libs[i]} $suffix ${spikes[i]}
done






2b. Normalize codon count files by rpm and rps

Inputs: 
${LibName}_Asite_${sizeRange}_CodonCoverageEach.txt
${LibName}_Asite_${sizeRange}_CodonCoverageEach_P.bedGraph
${LibName}_Asite_${sizeRange}_CodonCoverageEach_M.bedGraph
Outputs: 
${LibName}_Asite_${sizeRange}_CodonCoverageEach.rpm.txt
${LibName}_Asite_${sizeRange}_CodonCoverageEach.rps.txt
${LibName}_Asite_${sizeRange}_CodonCoverageEach_P.rpm.bedGraph
${LibName}_Asite_${sizeRange}_CodonCoverageEach_P.rps.bedGraph
${LibName}_Asite_${sizeRange}_CodonCoverageEach_M.rpm.bedGraph
${LibName}_Asite_${sizeRange}_CodonCoverageEach_M.rps.bedGraph

Exp="CHLtest"
Libs=("Hela_nodrug" "Hela_drug")
sizeRange="32to38"
suffix=Asite_${sizeRange}_CodonCoverageEach
spikes=(3000 3000) # If no spike ins here, just add this so it will run

for i in {0..1}
do
sbatch ./NormalizeCodonCount.sh ${Libs[i]} $suffix ${spikes[i]}
done

# Combine outputs
sizeRange="33to39"
Experiment="CHLtest"
Libs=("Hela_nodrug" "Hela_drug")
sizeRange="32to38"

versions=". .rpm. .rps."
for version in $versions
do
paste <(awk '{print $1"\t"$2"\t"$3}' ${Libs[0]}_Asite_${sizeRange}_CodonCoverageEach${version}txt) <(awk '{print $3}' ${Libs[1]}_Asite_${sizeRange}_CodonCoverageEach${version}txt) > ${Experiment}_Asite_${sizeRange}_CodonCoverageEach${version}txt
done






######### The above commands were written to run on bash terminal
######### Many of the below commands were written to run on zsh terminal, so there are minor differences





3. Call pauses in each sample using sliding window


3a. Calculate pause scores
### Input note ### 
### Input bedGraphs should be in a directory called ./bedGraphs/

Libs=("Hela_nodrug" "Hela_drug")
sizerange='32to38' 

version='rpm'
window=21
use='fold' 
# 'fold' 'SD' (use fold over mean or SD over mean to calc pause score?)

for lib in $Libs
do
Rscript ./Scripts/SlidingWindowPauseCaller.R $lib $sizerange $version $window $use
done



3b. Now call pauses from the sliding window pause scores
### Output note ### 
### Input bedGraphs should be in a directory called ./bedGraphs/

Experiment="CHLtest"
soi="Hela_drug" 
# If more than one, list as "Da,Db,Dc"
soiName="Drug" 
# short descriptive name of choice
sizeRange="32to38"
threshold=2 
use="fold"
length=21
# length is the number of nt on either side of the A site codon in the bed file, used to pull out sequence from fasta. Make it a factor of 3

Rscript ./Scripts/CallPausesFromSlidingWindows.R $Experiment $soi $soiName $sizeRange $threshold $use




3c. Pull out pause site sequence from fasta

infasta='./SeqFiles/hMitoGenome.fasta'
# samtools faidx $infasta
# sois="CHL LZD IND DMSO"
sois="Drug Nodrug"
cutoff='2' 
use='fold' # fold SD
suffix=slidingWindowPauses${use}_Thresh
for soi in $sois
do
bedtools getfasta -fo ${soi}_${suffix}${cutoff}_byCodon.fasta -s -name -fi $infasta -bed ${soi}_${suffix}${cutoff}_byCodon.bed
done





################ 4 #################
########## Motif logos #############
####################################

4. Get pause site logos by copy-pasting files into web apps:
# Translate here
https://www.bioinformatics.org/sms2/translate.html
# Get logos here
http://kplogo.wi.mit.edu/

############### end 4 ##############
####################################







################ 5 #################
####### Read/pause counts ##########
####################################

5. Get fraction of reads and/or pauses in each sample that have codon in each ribosome position (e.g. P-2, E(P-1), P, A(P+1), P+2)

# First need to add in strand info and go from bedgraph to bed)

Libs=("Ca" "Cb" "Cc" "Ce" "Cf" "Da" "Db" "De" "Df" "La" "Lb" "Lc" "Ld" "Le" "Lf")
sizerange='33to39' 

for lib in $Libs
do
awk 'NR>1 {print $1"\t"$2-21"\t"$3+21"\t0\t"$4"\t+"}' ${lib}_Asite_${sizerange}_CodonCoverageEach_P.rpm.bedGraph > ${lib}_Asite_${sizerange}_CodonCoverageEach_P.rpm.bed
awk 'NR>1 {print $1"\t"$2-21"\t"$3+21"\t0\t"$4"\t-"}' ${lib}_Asite_${sizerange}_CodonCoverageEach_M.rpm.bedGraph > ${lib}_Asite_${sizerange}_CodonCoverageEach_M.rpm.bed
cat ${lib}_Asite_${sizerange}_CodonCoverageEach_P.rpm.bed ${lib}_Asite_${sizerange}_CodonCoverageEach_M.rpm.bed > ${lib}_Asite_${sizerange}_CodonCoverageEach.rpm.bed
done



# Now get fasta sequence

module load gcc/9.2.0
module load bedtools/2.30.0
infasta='./SeqFiles/hMitoGenome.fasta'

# For reads
Libs="Ca Cb Cc Ce Cf Da Db De Df La Lb Lc Ld Le Lf"
sizerange='33to39' 
for lib in $Libs
do
bedtools getfasta -bedOut -s -name -fi $infasta -bed ${lib}_Asite_${sizerange}_CodonCoverageEach.rpm.bed > ${lib}_Asite_${sizerange}_CodonCoverageEach.rpm.fasta
done

# For pauses
Libs="DMSO CHL LZD"
for lib in $Libs
do
bedtools getfasta -bedOut -s -name -fi $infasta -bed ${lib}_slidingWindowPausesfold_Thresh2_byCodon.bed > ${lib}_slidingWindowPausesfold_Thresh2_byCodon.tab.fasta
done


# Now make barplots for each amino acid for each position

# This script will have to be opened and modified (paths, colors, samples, positions, etc) for each new experiment
script="./Scripts/TranslateANDcount_aas.r"

type="slidingWindowPauses"
Rscript $script $type

############### end 5 ##############
####################################







################ 6 #################
####### Codon occupancies ##########
####################################

6. Get codon occupancies
# These scripts will make bar plots of codon occupancy for each codon. 

mkdir CodonFrequency
cd CodonFrequency

# List this way because running on zsh terminal here
Libs=("Hela_nodrug" "Hela_drug")
sizerange="32to38"

norm="rpm"

site='PA'

for lib in $Libs
do
echo $lib
~/Desktop/Data/hMitoRP/PanAnalysis/ScriptsForReanalysis/CodonOccupancy.sh $lib $sizerange $site $norm
done



# Combine outputs
Exp="CHLtest"
Libs=("Hela_nodrug" "Hela_drug")
sizeRange='32to38' 

site="PA"
# Can be E P A EP PA EA

# # Start with 1 here instead of 0 for zsh 
paste <(awk '{print $1"\t"$2"\t"$3}' ${Libs[1]}_${sizeRange}_${site}siteCodonDensities.txt) <(awk '{print $2"\t"$3"\t"$4}' ${Libs[2]}_${sizeRange}_${site}siteCodonDensities.txt) > ${Exp}_${sizeRange}_${site}siteCodonDensities.txt


# open and *** modify *** then run R scripts
mkdir plots

# For a single site (e.g. A or P)
Experiment="CHLtest"
site="P"
soiName="Drug"
cutoff=4
#cutoff is how many times the codon must be present in the transcriptome to be considered
Rscript ./Scripts/CodonDensities_tRNAabundance.R $Experiment $site $soiName $cutoff

# For dicodon site (e.g. EP)
Experiment="CHLtest"
site="EP"
soiName="Drug"
cutoff=2
Rscript ./Scripts/diCodonDensities_tRNAabundance_foldchange.R $Experiment $site $soiName $cutoff

############### end 6 ##############
####################################







