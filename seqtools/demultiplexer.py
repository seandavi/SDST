#!/usr/bin/env python
import argparse
from itertools import izip
import gzip
from string import maketrans
import subprocess
import os
from multiprocessing import Pool
import collections


transtab = maketrans('ACGTN','TGCAN')

revcompdict = {}
def revcomp(sequence):
    """Reverse complement a string

    :param sequence: The DNA string (all caps)
    """
    if(sequence in revcompdict):
        return revcompdict[sequence]
    tmp=sequence[::-1].translate(transtab)
    revcompdict[sequence]=tmp
    return(tmp)


def fileOpen(fname,mode='r'):
    """Open a file, including gzip files

    :param fname: The filename to open.  Gzip files are distinguished by ending in '.gz'
    :param mode: File mode for opening [default 'r']
    :returns: An open file handle"""
    # gzip in python is REALLY slow, so use pipes instead.
    if(fname.endswith('.gz')):
        if(mode.startswith('r')):
            return subprocess.Popen(['gunzip -c %s' % fname],stdout=subprocess.PIPE,shell=True).stdout
        if(mode.startswith('w')):
            return subprocess.Popen(['gzip > %s' % fname],stdin=subprocess.PIPE,shell=True).stdin
    else:
        return open(fname,'r')



class FastqRecord(object):
    """Very simple fastq class containing header, sequence, line3, and quality as strings"""
    def __init__(self,header,sequence,line3,quality):
        """Initialize a FastqRecord

        :param header: header
        :param sequence: sequence
        :param line3: second header line from fastq record
        :param quality: string quality"""
        self.header=header
        self.sequence=sequence
        self.line3=line3
        self.quality=quality

    def __str__(self):
        """String represendation, suitable for output to a fastq file"""
        return "%s\n%s\n%s\n%s\n" % (self.header,self.sequence,self.line3,self.quality)
    

class Fastq(object):
    def __init__(self,fname):
        self.name = fname
        self.fh = fileOpen(fname)
        
    def _getNextRecord(self):
        x = []
        for i in range(0,4):
            nextLine=self.fh.next().strip()
            x.append(nextLine)
            if(nextLine==''):
                return None
        return(FastqRecord(x[0],x[1],x[2],x[3]))

    def __iter__(self):
        return self

    def next(self):
        tmp = self._getNextRecord()
        if(tmp is None):
            raise StopIteration
        else:
            return(tmp)

    def name(self):
        return self.fname

    def close(self):
        self.fh.close()
    
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
        

