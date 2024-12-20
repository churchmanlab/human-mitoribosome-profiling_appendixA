#!/bin/bash

#SBATCH -c 1
#SBATCH -t 0-0:05
#SBATCH -p short
#SBATCH --mem=100M
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mtcouvi@gmail.com



LibName=$1
suffix=$2 # Mito_mRNA.noDups.5p.plus Mito_mRNA.noDups.Asite_${sizeRange}_P
spike=$3


infile=${LibName}_${suffix}.txt
outfile=${LibName}_${suffix}.rpm.txt
outfilerps=${LibName}_${suffix}.rps.txt

# Get sum of all reads
sum=`awk 'BEGIN {OFS = "\t"} {sum+= $3; next;} END {print sum}' ${infile}`
# Normalize rpm
head -1 $infile > $outfile
awk 'FNR == 1 {next} BEGIN {OFS = "\t"} {print $1, $2, $3/s*1000000}' s=$sum $infile >> $outfile
# Normalize rps
head -1 $infile > $outfilerps
awk 'FNR == 1 {next} BEGIN {OFS = "\t"} {print $1, $2, $3/s*10000}' s=$spike $infile >> $outfilerps


infileP=${LibName}_${suffix}_P.bedGraph
infileM=${LibName}_${suffix}_M.bedGraph

outfileP=${LibName}_${suffix}_P.rpm.bedGraph
outfilerpsP=${LibName}_${suffix}_P.rps.bedGraph

outfileM=${LibName}_${suffix}_M.rpm.bedGraph
outfilerpsM=${LibName}_${suffix}_M.rps.bedGraph

# Get sum of all reads
# Get sum of all reads
sumP=`awk 'BEGIN {OFS = "\t"} {sum+= $4; next;} END {print sum}' ${infileP}`
sumM=`awk 'BEGIN {OFS = "\t"} {sum+= $4; next;} END {print sum}' ${infileM}`
sum=`expr $sumP + $sumM`

# Normalize rpm
head -1 $infileP > $outfileP
awk 'FNR == 1 {next} BEGIN {OFS = "\t"} {print $1, $2, $3, $4/s*1000000}' s=$sum $infileP >> $outfileP
head -1 $infileM > $outfileM
awk 'FNR == 1 {next} BEGIN {OFS = "\t"} {print $1, $2, $3, $4/s*1000000}' s=$sum $infileM >> $outfileM
# Normalize rps
head -1 $infileP > $outfilerpsP
awk 'FNR == 1 {next} BEGIN {OFS = "\t"} {print $1, $2, $3, $4/s*10000}' s=$spike $infileP >> $outfilerpsP
head -1 $infileM > $outfilerpsM
awk 'FNR == 1 {next} BEGIN {OFS = "\t"} {print $1, $2, $3, $4/s*10000}' s=$spike $infileM >> $outfilerpsM


