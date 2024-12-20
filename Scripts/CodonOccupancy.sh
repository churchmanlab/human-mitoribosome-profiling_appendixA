#!/bin/bash


### WRITTEN:		06/10/2020 by Mary Couvillion

### USE: 			This needs to be run from within 'CodonFrequency' directory
###					with 'bedGraphs' directory on the same level, containing (+)
###					and (-) bedGraph files

libName=$1
sizeRange=$2
site=$3
norm=$4

echo 'Starting get_gene_sequence_densities_MTC_hMitoRP.py'
python ../Scripts/get_gene_sequence_densities_MTC_hMitoRP.py -f ../SeqFiles/hMitoExonRanges.txt -g ../SeqFiles/hMitoGenome.fasta -l ${libName} -p Mito_mRNA.noDups.Asite_${sizeRange}_P.${norm}.bedGraph -m Mito_mRNA.noDups.Asite_${sizeRange}_M.${norm}.bedGraph -y $site
echo 'Done getting gene sequence densities'

echo 'Starting calculate_codon_density_MTC_hMitoRP.py'
python ../Scripts/calculate_codon_density_MTC_hMitoRP.py -f MitoGenes_seq.txt -s ${sizeRange} -y $site -l ${libName} -c ../SeqFiles/hMitoGeneticCode_BattersbyMyoblast_tRNAabundance.txt 
echo 'Done calculating codon densities'


# Add header (need all this extra stuff because of OSX, wouldn't if on O2)
sed -i -e '1i\'$'\n'$'codon'$'\t'${libName}'_aveDensity'$'\t'${libName}'_SD'$'\t''CodonFreq'$'\n' ${libName}_${sizeRange}_${site}siteCodonDensities.txt
