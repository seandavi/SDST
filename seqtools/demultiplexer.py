#!/usr/bin/env python
import argparse
from itertools import izip
import subprocess
import os
import collections
from seqtools.utils import revcomp
from seqtools.fastq import Fastq


def isIndexRevComp(indexfile,indexes,n=500000):
    """Determine if the indexes are reverse complemented or not
    
    :param indexfile: filename of the Fastq index file
    :param indexes: list or tuple of index strings
    :param n: integer number of reads to sample
    """
    ifile = Fastq(indexfile)
    indexreads = collections.defaultdict(int)
    for i in xrange(n):
        indexreads[ifile.next().sequence]+=1
    counts = {'normal':0,
              'revcomp':0}
    for k,v in indexreads.items():
        for i in indexes:
            if(k.startswith(i)):
                counts['normal']+=v
            if(k.startswith(revcomp(i))):
                counts['revcomp']+=v
    return(counts['revcomp']>counts['normal'])
    

def demultiplex(readfile,
                indexfile,
                indexes,
                readfile2=None,
                indexfile2=None,
                indexes2=None):

    # single readfile, single indexfile
    if(readfile2 is None) and (indexfile2 is None):
        rfile1 = Fastq(readfile)
        (rpath,rname) = os.path.split(readfile)
        ifile = Fastq(indexfile)
        indexRevComp = isIndexRevComp(indexfile,indexes)
        existingIndexes = []
        for i in indexes:
            ofname1 = os.path.join(rpath,i + "_" + rname)
            if(not os.path.exists(ofname1)):
                ofile1[i]=fileOpen(os.path.join(rpath,"tmp." + i + "_" + rname),'w')
            else:
                print ofname1," already exists, skipping"
                existingIndexes.append(i)
        for i in existingIndexes:
            indexes.remove(i)
        if(len(indexes)==0):
            exit(0)
        for (r1,i) in izip(rfile1,ifile):
            try:
                if indexRevComp:
                    i2 = revcomp(i.sequence[:indexlen])
                    ofile1[i2].write(str(r1))
                else:
                    i2 = i.sequence[:indexlen]
                    ofile1[i2].write(str(r1))
            except KeyError:
                pass
        rfile1.close()
        ifile.close()
        for ofile in ofile1.values():
            ofile.close()
        for i in indexes:
            os.rename(os.path.join(rpath,'tmp.' + i + "_" + rname),
                      os.path.join(rpath,i + "_" + rname))
            
    # two readfiles, single indexfile
    if(readfile2 is not None) and (indexfile2 is None):
        rfile1 = Fastq(readfile)
        rfile2 = Fastq(readfile2)
        (rpath,rname) = os.path.split(readfile)
        (rpath2,rname2) = os.path.split(readfile2)
        ifile = Fastq(indexfile)
        indexRevComp = isIndexRevComp(indexfile,indexes)
        ofile1 = {}
        ofile2 = {}
        existingIndexes = []
        for i in indexes:
            ofname1 = os.path.join(rpath,i + "_" + rname)
            ofname2 = os.path.join(rpath2,i + "_" + rname2)
            if(os.path.exists(ofname1) and os.path.exists(ofname2)):
                print ofname1,ofname2, " already exist, skipping"
                existingIndexes.append(i)
            else:
                ofile1[i]=fileOpen(os.path.join(rpath,"tmp." + i + "_" + rname),'w')
                ofile2[i]=fileOpen(os.path.join(rpath2,"tmp." + i + "_" + rname2),'w')
        for i in existingIndexes:
            indexes.remove(i)
        if(len(indexes)==0):
            exit(0)
        indexlen = len(indexes[0])
        for (r1,r2,i) in izip(rfile1,rfile2,ifile):
            try:
                if indexRevComp:
                    i2 = revcomp(i.sequence[:indexlen])
                    ofile1[i2].write(str(r1))
                    ofile2[i2].write(str(r2))
                else:
                    i2 = i.sequence[:indexlen]
                    ofile1[i2].write(str(r1))
                    ofile2[i2].write(str(r2))                    
            except KeyError:
                pass
        rfile1.close()
        rfile2.close()
        ifile.close()
        for ofile in ofile1.values():
            ofile.close()
        for ofile in ofile2.values():
            ofile.close()
        for i in indexes:
            print os.path.join(rpath,'tmp.' + i + "_" + rname),                      os.path.join(rpath,i + "_"+rname)
            os.rename(os.path.join(rpath,'tmp.' + i + "_" + rname),
                      os.path.join(rpath,i + "_"+rname))
            os.rename(os.path.join(rpath2,'tmp.' + i +"_"+ rname2),
                      os.path.join(rpath2,i +"_"+ rname2))

def main():
    demultiplex(readfile='/data/CCRBioinfo/fastq/3_1_AD0P2BACXX.319_BUSTARD-2012-03-17.fq.gz',
                indexfile='/data/CCRBioinfo/fastq/3_2_AD0P2BACXX.319_BUSTARD-2012-03-17.fq.gz',
                indexes=('TGGTCA','CACTGT','ATTGGC'),
                readfile2='/data/CCRBioinfo/fastq/3_3_AD0P2BACXX.319_BUSTARD-2012-03-17.fq.gz')
    
    exit(0)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-1','--readFile1',
                        help="read1 filename")
    parser.add_argument('-2','--readFile2',
                        help="read2 filename")
    parser.add_argument('-i','--indexFile',
                        help="index Filename")
    parser.add_argument('-x','--index',type=str,action='append',
                        help="The indexes, one per index")

    opts = parser.parse_args()


    if(opts.readFile2):
        demultiplex(readfile = opts.readFile1,
                    readfile2 = opts.readFile2,
                    indexfile = opts.indexFile,
                    indexes = opts.index)
    else:
        demultiplex(readfile = opts.readFile1,
                    indexfile = opts.indexFile,
                    indexes = opts.index)
        

