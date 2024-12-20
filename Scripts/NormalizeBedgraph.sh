#!/bin/bash

#SBATCH -c 1
#SBATCH -t 0-0:05
#SBATCH -p short
#SBATCH --mem=1G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mtcouvi@gmail.com



LibName=$1
suffix=$2 # Mito_mRNA.noDups.5p.plus Mito_mRNA.noDups.Asite_${sizeRange}_P
spike=$3


# Plus
infileP1=${LibName}_${suffix}_P.bedGraph
outfileP1=${LibName}_${suffix}_P.rpm.bedGraph
outfilerpsP1=${LibName}_${suffix}_P.rps.bedGraph
# Minus
infileM1=${LibName}_${suffix}_M.bedGraph
outfileM1=${LibName}_${suffix}_M.rpm.bedGraph
outfilerpsM1=${LibName}_${suffix}_M.rps.bedGraph


# Plus
infileP2=${LibName}_${suffix}_PAll.txt
outfileP2=${LibName}_${suffix}_PAll.rpm.bedGraph
outfilerpsP2=${LibName}_${suffix}_PAll.rps.bedGraph
# Minus
infileM2=${LibName}_${suffix}_MAll.txt
outfileM2=${LibName}_${suffix}_MAll.rpm.bedGraph
outfilerpsM2=${LibName}_${suffix}_PAll.rps.bedGraph


# Get sum of all reads
sumP=`awk 'BEGIN {OFS = "\t"} {sum+= $4*($3-$2); next;} END {print sum}' ${infileP1}`
sumM=`awk 'BEGIN {OFS = "\t"} {sum+= $4*($3-$2); next;} END {print sum}' ${infileM1}`
sum=`expr $sumP + $sumM`
# Normalize to rpm
awk 'BEGIN {OFS = "\t"} {print $1, $2, $3, $4/s*1000000}' s=$sum $infileP1 > $outfileP1
awk 'BEGIN {OFS = "\t"} {print $1, $2, $3, $4/s*1000000}' s=$sum $infileM1 > $outfileM1
# Normalize file with all positions filled
awk 'BEGIN {OFS = "\t"} {print $1, $2, $3, $4/s*1000000}' s=$sum $infileP2 > $outfileP2
awk 'BEGIN {OFS = "\t"} {print $1, $2, $3, $4/s*1000000}' s=$sum $infileM2 > $outfileM2

# Normalize to rps
awk 'BEGIN {OFS = "\t"} {print $1, $2, $3, $4/s*10000}' s=$spike $infileP1 > $outfilerpsP1
awk 'BEGIN {OFS = "\t"} {print $1, $2, $3, $4/s*10000}' s=$spike $infileM1 > $outfilerpsM1
# Normalize file with all positions filled
awk 'BEGIN {OFS = "\t"} {print $1, $2, $3, $4/s*10000}' s=$spike $infileP2 > $outfilerpsP2
awk 'BEGIN {OFS = "\t"} {print $1, $2, $3, $4/s*10000}' s=$spike $infileM2 > $outfilerpsM2
