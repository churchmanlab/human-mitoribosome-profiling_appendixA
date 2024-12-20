# Count reads in each frame for exons


#On command line write python CodonCoverage.py <plus bedGraph> <minus bedGraph> <LibName>
import sys
infilename_plus = sys.argv[1]
infilename_minus = sys.argv[2]
LibName = sys.argv[3]
CodonFile = open('./SeqFiles/cutCodons_hMito_noOverlap.bed','r') 
AllCodonFile = open('./SeqFiles/cutCodons_hMito_5pgene.bed','r') 


#Open the bedGraph file to read
InputFilePlus = open(infilename_plus, 'r')
InputFileMinus = open(infilename_minus, 'r')

#Make dictionary for positions and read counts
PlusCounts = {}
MinusCounts = {}

for line in InputFilePlus:
    if line[0:2] != 'tr':
        clist = line.replace('\n','').split('\t')
        PlusCounts[clist[2]] = [clist[3]]
for line in InputFileMinus:
    if line[0:2] != 'tr':
        clist = line.replace('\n','').split('\t')
        MinusCounts[clist[2]] = [clist[3]]
# InputFilePlus.seek(0)
# InputFileMinus.seek(0)

# Go through codon file and make list entry for counts at each codon
PlusCodonCounts = []
MinusCodonCounts = []
genes = []
starts = []
n=0

# First use noOverlap file to get overall coverage values
for line in CodonFile:
    n = n+1
    cols = line.replace('\n','').split('\t')
    genes.append(cols[4])
    starts.append(cols[1])
    POSs = range(int(cols[1]),int(cols[2])+1,1)
    POSs = [str(i) for i in POSs]
    if cols[3] == '+' and int(cols[4].split('.')[-1]) > 6:
        codonCount = float(PlusCounts[POSs[0]][0]) + float(PlusCounts[POSs[1]][0]) + float(PlusCounts[POSs[2]][0])
        PlusCodonCounts.append(codonCount)
    elif cols[3] == '-' and int(cols[4].split('.')[-1]) > 6:
        codonCount = float(MinusCounts[POSs[0]][0]) + float(MinusCounts[POSs[1]][0]) + float(MinusCounts[POSs[2]][0])
        MinusCodonCounts.append(codonCount)

# Get total number of codons and number with >0
CodonCounts = MinusCodonCounts + PlusCodonCounts
Tot = len(CodonCounts)
Covered = Tot - CodonCounts.count(0)
PerCovered = round(float(Covered)/float(Tot) * 100, 3)

def mean(lst):
    return sum(lst) / len(lst)

MeanCoverage = mean(CodonCounts)


#Create a file to write
OutputFile = open((LibName + '_CodonCoverage_ignore1st6.txt'), 'w')

OutputFile.write('Percent codons covered' + '\t' + str(PerCovered) + '\n' + 'Average codon coverage' + '\t' + str(MeanCoverage))







# Now use file with all the codons included to get individual coverage
PlusCodonCounts = []
MinusCodonCounts = []
Plusgenes = []
Plusstarts = []
Minusgenes = []
Minusstarts = []
p=0
m=0

for line in AllCodonFile:
    cols = line.replace('\n','').split('\t')
    if cols[3] == '+':
        p=p+1
        Plusgenes.append(cols[4])
        Plusstarts.append(cols[1])
        POSs = range(int(cols[1]),int(cols[2])+1,1)
        POSs = [str(i) for i in POSs]
        codonCount = float(PlusCounts[POSs[0]][0]) + float(PlusCounts[POSs[1]][0]) + float(PlusCounts[POSs[2]][0])
        PlusCodonCounts.append(codonCount)
    elif cols[3] == '-':
        m=m+1
        Minusgenes.append(cols[4])
        Minusstarts.append(cols[1])
        POSs = range(int(cols[1]),int(cols[2])+1,1)
        POSs = [str(i) for i in POSs]
        codonCount = float(MinusCounts[POSs[0]][0]) + float(MinusCounts[POSs[1]][0]) + float(MinusCounts[POSs[2]][0])
        MinusCodonCounts.append(codonCount)



OutputFile2 = open((LibName + '_CodonCoverageEach.txt'), 'w')
OutputFile2.write('Pos\tntStart\tCounts\n')
for i in range(0,p):
    print(i)
    OutputFile2.write(Plusgenes[i] + '\t' + str(Plusstarts[i]) + '\t' +  str(PlusCodonCounts[i]) + '\n')
for i in range(0,m):
    print(i)
    OutputFile2.write(Minusgenes[i] + '\t' + str(Minusstarts[i]) + '\t' +  str(MinusCodonCounts[i]) + '\n')

# And one to load on IGV
OutputFile3P = open((LibName + '_CodonCoverageEach_P.bedGraph'), 'w')
OutputFile3M = open((LibName + '_CodonCoverageEach_M.bedGraph'), 'w')
OutputFile3P.write('track type=bedGraph\n')
for i in range(0,p):
    OutputFile3P.write('chrM' + '\t' + str(int(Plusstarts[i])-1) + '\t' + str(int(Plusstarts[i])+2) + '\t' +  str(PlusCodonCounts[i]) + '\n')
OutputFile3M.write('track type=bedGraph\n')
for i in range(0,m):
    OutputFile3M.write('chrM' + '\t' + str(int(Minusstarts[i])-1) + '\t' + str(int(Minusstarts[i])+2) + '\t' +  str(MinusCodonCounts[i]) + '\n')




InputFilePlus.close()
InputFileMinus.close()
OutputFile.close()

