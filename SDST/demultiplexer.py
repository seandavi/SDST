import argparse

import subprocess
import os
import collections
import pylev

from SDST.utils import revcomp,fileOpen
from SDST.fastq import Fastq

def isIndexRevComp(indexfile,indexes,n=500000):
    """Determine if the indexes are reverse complemented or not
    
    :param indexfile: filename of the Fastq index file
    :param indexes: list or tuple of index strings
    :param n: integer number of reads to sample
    """
    print("HERE")
    ifile = Fastq(indexfile)
    ilength=len(indexes[0])
    print(ilength)
    indexreads = collections.defaultdict(int)
    for i in range(n):
        indexreads[ifile.next().sequence[:ilength]]+=1
    counts = {'normal':0,
              'revcomp':0}
    for k,v in list(indexreads.items()):
        print(k,v)
        for i in indexes:
            if(pylev.levenshtein(k,i)<=1):
                counts['normal']+=v
                continue
            if(pylev.levenshtein(k,revcomp(i))<=1):
                counts['revcomp']+=v
    if(counts['revcomp']>counts['normal']):
        print('using revcomp')
    else:
        print('NOT revcomp')
        
    return(counts['revcomp']>counts['normal'])
    

def demultiplex(readfile,
                indexfile,
                indexes,
                readfile2=None,
                indexfile2=None):
    """Demultiplex from separate FASTQ files.

    All FASTQ files can be gzipped (with suffix .gz).

    :param readfile: The filename of the first fastq file
    :param indexfile: The filename of the first index fastq file
    :param indexes: An iterable of indexes.  If dual-barcoding is used, the indexes should be comma-separated strings, one string for each barcode pair.
    :param indexfile2: The filename of the second index fastq file.  If this parameter is included, then the indexes parameter should be a set of comma-separated pairs of indexes.  
    :param readfile2: The filename of the second fastq file [optional]
    
    """

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
                ofile1[i]=fileOpen(os.path.join(rpath,i + "_" + rname),'w')
            else:
                print(ofname1," already exists, skipping")
                existingIndexes.append(i)
        for i in existingIndexes:
            indexes.remove(i)
        if(len(indexes)==0):
            exit(0)
        for (r1,i) in zip(rfile1,ifile):
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
        for ofile in list(ofile1.values()):
            ofile.close()
        ## for i in indexes:
        ##     os.rename(os.path.join(rpath,'tmp.' + i + "_" + rname),
        ##               os.path.join(rpath,i + "_" + rname))
            
    # two readfiles, single indexfile
    if(readfile2 is not None) and (indexfile2 is None):
        print("here1")
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
                print(ofname1,ofname2, " already exist, skipping")
                existingIndexes.append(i)
            else:
                ofile1[i]=fileOpen(os.path.join(rpath,i + "_" + rname),'w')
                ofile2[i]=fileOpen(os.path.join(rpath2,i + "_" + rname2),'w')
        for i in existingIndexes:
            indexes.remove(i)
        if(len(indexes)==0):
            exit(0)
        indexlen = len(indexes[0])
        for (r1,r2,i) in zip(rfile1,rfile2,ifile):
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
        for ofile in list(ofile1.values()):
            ofile.close()
        for ofile in list(ofile2.values()):
            ofile.close()
        ## for i in indexes:
        ##     print os.path.join(rpath,'tmp.' + i + "_" + rname),os.path.join(rpath,i + "_"+rname)
        ##     os.rename(os.path.join(rpath,'tmp.' + i + "_" + rname),
        ##               os.path.join(rpath,i + "_"+rname))
        ##     os.rename(os.path.join(rpath2,'tmp.' + i +"_"+ rname2),
        ##               os.path.join(rpath2,i +"_"+ rname2))

    # two readfiles, two indexfiles
    if(readfile2 is not None) and (indexfile2 is not None):
        rfile1 = Fastq(readfile)
        rfile2 = Fastq(readfile2)
        (rpath,rname) = os.path.split(readfile)
        (rpath2,rname2) = os.path.split(readfile2)
        ifile = Fastq(indexfile)
        ifile2 = Fastq(indexfile2)
        indexes = [tuple(x.split(',')) for x in indexes]
        indexRevComp = isIndexRevComp(indexfile,[i[0] for i in indexes])
        ofile1 = {}
        ofile2 = {}
        existingIndexes = []
        for j in indexes:
            i = ''.join(j)
            ofname1 = os.path.join(rpath,i + "_" + rname)
            ofname2 = os.path.join(rpath2,i + "_" + rname2)
            if(os.path.exists(ofname1) and os.path.exists(ofname2)):
                print(ofname1,ofname2, " already exist, skipping")
                existingIndexes.append(i)
            else:
                ofile1[i]=fileOpen(ofname1,'w')
                ofile2[i]=fileOpen(ofname2,'w')
        for i in existingIndexes:
            indexes.remove(i)
        if(len(indexes)==0):
            exit(0)
        indexlen = len(indexes[0][0])
        for (r1,r2,i,i2) in zip(rfile1,rfile2,ifile,ifile2):
            try:
                if indexRevComp:
                    ir = revcomp(i.sequence[:indexlen])
                    ir2 = revcomp(i2.sequence[:indexlen])
                    istr = ir+ir2
                    ofile1[istr].write(str(r1))
                    ofile2[istr].write(str(r2))
                else:
                    ir = i.sequence[:indexlen]
                    ir2 = i2.sequence[:indexlen]
                    istr = ir+ir2
                    ofile1[istr].write(str(r1))
                    ofile2[istr].write(str(r2))
            except KeyError:
                pass
        rfile1.close()
        rfile2.close()
        ifile.close()
        ifile2.close()
        for ofile in list(ofile1.values()):
            ofile.close()
        for ofile in list(ofile2.values()):
            ofile.close()
        ## for i in indexes:
        ##     ofname1 = os.path.join(rpath,''.join(i) + "_" + rname)
        ##     ofname2 = os.path.join(rpath2,''.join(i) + "_" + rname2)
        ##     os.rename(os.path.join(rpath,'tmp.' + ofname1),
        ##               os.path.join(rpath,ofname1))
        ##     os.rename(os.path.join(rpath2,'tmp.'+ofname2),
        ##               os.path.join(rpath2,ofname2))


