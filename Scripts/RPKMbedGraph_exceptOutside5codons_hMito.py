HELP_STRING = """
RPKM.py

Load a (pre-normalized) plus and minus bedGraph file to get RPKM values

Author: Mary Couvillion
Date: 3/3/15
Modified 10/2016 for human genome


     -h     print this help message
     -p     plus strand bedGraph (required)
     -m     minus strand bedGraph (required)
     -o     output file name (required)
     
command line, e.g.: python ~/Desktop/Data/Python\ scripts/RPKMfrom_bedGraph.py -p CP11_AsitePlusRPM_Norm.bedGraph -m CP11_AsiteMinusRPM_Norm.bedGraph -o CP11_RPKM_new.txt

for running on Orchestra:  bsub -q short -W 12:00 -o logs/CP17_RPKM.log -e logs/CP17_RPKM.err python /home/mc348/RPKMbedGraph_exceptOutside5codons.py -p CP17_Nuc_mRNA/CP17_AsitePlusRPM_Norm.bedGraph -m CP17_Nuc_mRNA/CP17_AsiteMinusRPM_Norm.bedGraph -o CP17_RPKM_inside.txt

"""
import sys
import itertools
from getopt import getopt


def main(argv=None):
    if argv is None:
        argv = sys.argv

    plusFile = ""
    minusFile = ""
    outputFile = ""

    try:
        optlist, args = getopt(argv[1:], "hp:m:o:")
    except:
        print ""
        print HELP_STRING
        sys.exit(1)
       
    if len(optlist) == 0:
        print ""
        print HELP_STRING
        sys.exit(1)
       
    for (opt, opt_arg) in optlist:
        if opt == '-h':
            print ""
            print HELP_STRING
            sys.exit(1)
        elif opt == '-p':
            plusFile = opt_arg
        elif opt == '-m':
            minusFile = opt_arg
        elif opt == '-o':
            outputFile = opt_arg
    if plusFile == '' or minusFile == '' or outputFile == "":
        print HELP_STRING
        sys.exit(1)

    PlusFile = open(plusFile,'r')
    MinusFile = open(minusFile, 'r')
    InputFile = open('/home/mc348/hMitoExonRanges_noOverlap.txt', 'r') #/Users/Mary/Desktop/Data/hMitoRP/SequenceFiles/hMitoExonRanges.txt
    
    outFile = open(outputFile, 'w')

    outFile.write('%s\t%s\t%s\t%s\t%s\n' % ('gene_name', 'strand', 'locus', 'RPK', 'Length'))

    for i,line in enumerate(InputFile):
        if i%200 == 0:
            print i, line

        clist = line.replace('\n','').split('\t')
        exon_ranges = clist[4].split(';')
        exon_number = len(exon_ranges)
        exon_ranges = [i.split('-') for i in exon_ranges]
        exon_ranges = [map(lambda x: int(x), i) for i in exon_ranges]
        positions_list = []
        for i in range(0,int(exon_number)):
            positions_list.append(range(exon_ranges[i][0], exon_ranges[i][1]+1))
        positions_all = list(itertools.chain.from_iterable(positions_list)) #flatten nested list
        positions=positions_all[15:-15]

        if clist[2] == '+':
            reads = []
            pos = -1
            lines_list = PlusFile.readlines()
            for row in lines_list:
                pos = pos+1
                if row[0] == 't':
                    continue
                Bed_clist = row.split()
                if int(Bed_clist[2]) in positions:
                    reads.append((float(Bed_clist[3]))*(float(Bed_clist[2])-float(Bed_clist[1])))

            RPM = sum(reads)
            RPK = RPM/((float(clist[3])-30)/1000) #adjust size for the 5 amino acids not counted at ends
            locus = clist[0] + ':' + str(positions_all[0]) + '-' + str(positions_all[-1])
            outFile.write('%s\t%s\t%s\t%s\t%s\n' % (clist[1], clist[2], locus,  RPK, clist[3]))
            PlusFile.seek(0)

        elif clist[2] == '-':
            reads = []
            pos=-1
            lines_list = MinusFile.readlines()
            for row in lines_list:
                pos = pos+1
                if row[0] == 't':
                    continue
                Bed_clist = row.split()
                if int(Bed_clist[2]) in positions:
                    reads.append((float(Bed_clist[3]))*(float(Bed_clist[2])-float(Bed_clist[1])))
            RPM = sum(reads)
            RPK = RPM/((float(clist[3])-30)/1000) #adjust size for the 5 amino acids not counted at ends
            locus = clist[0] + ':' + str(positions_all[0]) + '-' + str(positions_all[-1])
            outFile.write('%s\t%s\t%s\t%s\t%s\n' % (clist[1], clist[2], locus,  RPK, clist[3]))
            MinusFile.seek(0)

##############################################
if __name__ == "__main__":
    sys.exit(main())
