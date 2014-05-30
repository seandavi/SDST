#!/usr/bin/env python
import argparse
import shutil
import subprocess
import os

parser = argparse.ArgumentParser()
parser.add_argument('tumorbam',help='full path to the Tumor BAM file')
parser.add_argument('normalbam',help='full path to the Normal BAM file')
parser.add_argument('twobitfile',help='full path to 2bit file for reference genome')
parser.add_argument('fastafile',help='full path to fasta file for reference genome')
parser.add_argument('resultfile',help='The full path to the result file')
parser.add_argument('-w','--windowsize',default=10000000,type=int,
                    help='The windowsize for breaking up work')

args = parser.parse_args()

tumorbam = os.path.basename(args.tumorbam)
normalbam = os.path.basename(args.normalbam)
subprocess.call('clearscratch',shell=True)
shutil.copy(args.tumorbam,'/scratch/'+tumorbam)
shutil.copy(args.normalbam,'/scratch/'+normalbam)
shutil.copy(args.tumorbam+".bai",'/scratch/'+tumorbam+".bai")
shutil.copy(args.normalbam+".bai",'/scratch/'+normalbam+".bai")
subprocess.call('module load CREST; gfServer start localhost 50000 {0} -canStop -log=blatServer.log &'.format(args.twobitfile),shell=True)

def makeRegions(faifile,windowsize):
    """Given a .fai file, make windows of size windowsize as a generator, returning a tuple of (chrom,start,end)"""
    regions=[]
    with open(faifile,'r') as f:
        for row in f:
            srow = row.strip().split("\t")
            for start in range(1,int(srow[1]),windowsize):
                if((start+windowsize)>int(srow[1])):
                   regions.append((srow[0],start,int(srow[1])))
                else:
                   regions.append((srow[0],start,start+windowsize-1))
    return(regions)


def region2string(region):
    """Convert three-item tuple into chr1:1234-5678 format"""
    return('{0}:{1}-{2}'.format(region[0],region[1],region[2]))

def region2dotstring(region):
    """Convert three-item tuple into chr1.1234.5678 format used by CREST"""
    return('{0}.{1}.{2}'.format(region[0],region[1],region[2]))

def runCrest(region):
    """Run extractSClip and CREST on a region tuple (chrom,start,end)"""
    print(region)
    strregion=region2string(region)
    dotregion=region2dotstring(region)
    shellcmd = 'module load CREST; mkdir -p /scratch/{5}; cd /scratch; extractSClip.pl -r {1} -i {0} --ref_genome {3}; extractSClip.pl -r {1} -i {2} --ref_genome {3}; CREST.pl -f {0}.{5}.cover -d {0} -g {2} -r {1} -o /scratch/{5} --ref_genome {3} -t {4} --blatserver localhost --sensitive --max_rep_cover 750'.format(tumorbam,strregion,normalbam,args.fastafile,args.twobitfile,dotregion)
    print(shellcmd)
    subprocess.call(shellcmd,shell=True)


import multiprocessing
pool = multiprocessing.Pool(multiprocessing.cpu_count())
# chromosomes hard-coded for human and excludes random, etc.
pool.map(runCrest,makeRegions(args.fastafile+".fai",args.windowsize))
resfile = open(args.resultfile,'w')
for i in makeRegions(args.fastafile+".fai",args.windowsize):
    with open('/scratch/{0}/{1}.predSV.txt'.format(region2dotstring(region),tumorbam)) as f:
        for line in f:
            resfile.write(line)
resfile.close()
subprocess.call('gfServer stop localhost 50000',shell=True)






