plHELP_STRING = """
SAM2centerBED.py

Author: Mary Couvillion
Date: December 19, 2016 / A site transformation changed 9/25/18 (MC) / A site transformation changed again 6/2019
Date: April 29, 2020 Modified to include soft-clipped reads
Date: July 11, 2024 Modified to take lists of sizes and offsets
Here we take alignments for reads and convert them to A site positions. Then the alignments are put back into a BED format.
To interpret the SAM file: http://samtools.sourceforge.net/SAMv1.pdf and http://picard.sourceforge.net/explain-flags.html

NOTE this will not work for reads that span exon-exon junctions


     -h     print this help message
     -i     file input (required)
     -s    sizes e.g. "30 31 32"
     -o		offsets e.g. "13 14 15"
"""
import sys
from getopt import getopt


def main(argv=None):
    if argv is None:
        argv = sys.argv

    inputFile = ""
    outputFilePlus = ""
    outputFileMinus = ""
    sizes = ""
    offsets = ""
    
    try:
        optlist, args = getopt(argv[1:], "hi:s:o:")
    except:
        print ""
        print HELP_STRING
        sys.exit(1)
       
    if len(optlist) == 0:
        print ""
        print HELP_STRING
        sys.exit(1)
       
    for (opt, opt_arg) in optlist:
        #print opt
        #print opt_arg
        if opt == '-h':
            print ""
            print HELP_STRING
            sys.exit(1)
        elif opt == '-i':
            inputFile = opt_arg
        elif opt == '-s':
            sizes = opt_arg
        elif opt == '-o':
            offsets = opt_arg

    if inputFile == "" :
        print 'Here is the problem'
        sys.exit(1)
                
    
    def getAsitePlus(fivePrime,readSize, sizeVec, offsetVec):
        if read_size in sizeVec:
            threeP = fivePrime + readSize - 1
            ind=sizeVec.index(read_size)
            offset=offsetVec[ind]
            Asite = threeP - offset
        else:
            Asite = 'OutOfRange'
        return Asite

    def getAsiteMinus(threePrime,readSize,sizeVec, offsetVec):
        if read_size in sizeVec:
            ind=sizeVec.index(read_size)
            offset=offsetVec[ind]
            Asite = threePrime + offset
        else:
            Asite = 'OutOfRange'
        return Asite


    sizeRange=str(sizes[0:2]) + 'to' + str(sizes[-2:])
    InFile = file(inputFile,'r')
    outFilePlus = open(inputFile.replace('_for'+sizeRange+'.sam', '.Asite_'+sizeRange+'_P.bed'), 'w')
    outFileMinus = open(inputFile.replace('_for'+sizeRange+'.sam', '.Asite_'+sizeRange+'_M.bed'), 'w')
    flag = 0
    sizeVec=sizes.split()
    sizeVec=[int(i) for i in sizeVec]
    print sizeVec
    offsetVec=offsets.split()
    offsetVec=[int(i) for i in offsetVec]
    print offsetVec
    # iterate through and process all lines in input file
    for i,line in enumerate(open(inputFile)):
        if i%100000 == 0:
            print i
        
        clist = line.replace('\n','').split('\t')
        # only process line if the FLAG is 0 (plus strand primary) or 256 (plus strand not primary)
        if clist[0] != '@':
    
            if clist[1] == '0' or clist[1] == '256':

                #pull out the 5' coordinate of the sequening read
                fiveP = int(clist[3])
                # pull out the CIGAR string
                cigar = clist[5]
                # for now ignore reads with indels and split alignments
                if 'D' not in cigar and 'I' not in cigar and 'N' not in cigar:
                    if 'S' not in cigar:
                        # parse CIGAR string
                        read_size = int(cigar.split('M')[0])
                        Asite = getAsitePlus(fiveP,read_size,sizeVec, offsetVec)
                    # if there is soft clipping only at 5' end
                    elif cigar[-1] != 'S' and (cigar[1] == 'S' or cigar[2] == 'S'):
                        # Split by S then split by M
                        read_size = int(cigar.split('S')[1].split('M')[0])
                        Asite = getAsitePlus(fiveP,read_size,sizeVec, offsetVec)
                    # if there is soft clipping at the 3' end (doesn't matter whether 5' end is soft-clipped or not because cigar.split('S')[-2] will give the M part of the string either way
                    elif cigar[-1] == 'S':
                        tailLength = int(cigar.split('M')[-1].split('S')[0])
                        seq = clist[9]
                        tail = seq[-tailLength:]
                        count = 0
                        for letter in tail:
                            if letter != 'A':
                                count = count + 1
                            elif letter == 'A':
                                count = count
                        # If the tail is all As, add that length 
                        if count == 0:
                            read_size = int(cigar.split('S')[-2].split('M')[0]) + tailLength
                            Asite = getAsitePlus(fiveP,read_size,sizeVec, offsetVec)
                        # If not, ignore tail
                        elif count != 0 :
                            read_size = int(cigar.split('S')[-2].split('M')[0])
                            Asite = getAsitePlus(fiveP,read_size,sizeVec, offsetVec)                            
                    else:
                        continue
                    if Asite != 'OutOfRange':
                        start = Asite - 1
                        # Now we output the adjusted read in BED format
                        outFilePlus.write(clist[2]+'\t%s\t%s\n' % (start, Asite))
                    else:
                        flag = flag + 1

                else:
                    continue

            
            
            elif clist[1] == '16' or clist[1] == '272':

                #now first position listed is 3' end, because it's the reverse complement
                threeP = int(clist[3])
                # pull out the CIGAR string
                cigar = clist[5]
                # parse CIGAR string
                if 'D' not in cigar and 'I' not in cigar and 'N' not in cigar:
                    if 'S' not in cigar:
                        # parse CIGAR string
                        read_size = int(cigar.split('M')[0])
                        Asite = getAsiteMinus(threeP,read_size,sizeVec, offsetVec)
                    # if there is soft clipping only at 5' end (now 5' end is the right side of cigar string)
                    elif cigar[1] != 'S' and cigar[2] != 'S' and cigar[-1] == 'S':
                        read_size = int(cigar.split('M')[0])
                        Asite = getAsiteMinus(threeP,read_size,sizeVec, offsetVec)
                    # if there is soft clipping at the 3' end
                    elif cigar[1] == 'S' or cigar[2] == 'S':
                        tailLength = int(cigar.split('S')[0])
                        seq = clist[9]
                        tail = seq[0:tailLength]
                        count = 0
                        for letter in tail:
                            if letter != 'T':
                                count = count + 1
                            elif letter == 'T':
                                count = count
                        # If the tail is all Ts, add that length 
                        if count == 0:
                            read_size = int(cigar.split('S')[1].split('M')[0]) + tailLength
                            # For minus strand have to adjust the start coordinate
                            threePtail = threeP - tailLength
                            Asite = getAsiteMinus(threePtail,read_size,sizeVec, offsetVec)
                        # If not, ignore tail
                        elif count != 0 :
                            read_size = int(cigar.split('S')[1].split('M')[0])
                            Asite = getAsiteMinus(threeP,read_size,sizeVec, offsetVec)
                    else:
                        continue
                    
                    if Asite != 'OutOfRange':
                        start = Asite - 1
                        # Now we output the adjusted read in BED format
                        outFileMinus.write(clist[2]+'\t%s\t%s\n' % (start, Asite))
                    else:
                        flag = flag + 1

                else:
                    continue
    print flag, 'records were out of size range'


##############################################
if __name__ == "__main__":
    sys.exit(main())
