#!/usr/bin/env python
import argparse

def convertChrom(gasvChrom):
    """Convert numeric GASV coordinates to UCSC chromosomes"""
    chrom = 'chr'+gasvChrom
    if(chrom=='chr23'): return 'chrX'
    if(chrom=='chr24'): return 'chrY'
    return(chrom)

def gasvLine2bedpeLine(gasvline):
    sgasv = gasvline.strip().split('\t')
    chrom1 = convertChrom(sgasv[1])
    chrom2 = convertChrom(sgasv[3])
    bedpe = tuple([chrom1] + sgasv[2].split(',') + [chrom2] + sgasv[4].split(',') + [sgasv[0],sgasv[5],'.','.',sgasv[6],sgasv[7]])
    return(bedpe)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Convert GASV output to BEDPE format')
    parser.add_argument('fname')
    parser.add_argument('-r','--readPairs',type=int,
                        help="Minumum number of read pairs to pass through the conversion process",default=0)
    parser.add_argument('-f','--noOverlappingEnds',action='store_false',default=True)
    parser.add_argument('-s','--slop',type=int,
                        help="Minimum span between end and start of intrachromosomal events",default=1000)
    args = parser.parse_args()
    with open(args.fname,'r') as f:
        f.next()
        for line in f:
            bedpe = gasvLine2bedpeLine(line)
            if(int(bedpe[7])>args.readPairs):
                if(args.noOverlappingEnds):
                    if((bedpe[0]==bedpe[3]) and (int(bedpe[2])+args.slop>int(bedpe[4]))):
                       continue
                print("\t".join(bedpe))
