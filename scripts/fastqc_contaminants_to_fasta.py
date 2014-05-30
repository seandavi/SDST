#!/usr/bin/env python
# The scythe trimmer uses a fasta file of contaminants
# and the fastqc contaminants file in:
#
#  $FASTQC_ROOT/Contaminants/contaminant_list.txt
#
# is the best list I have found.
#
# This script takes as input the contaminant_list file
# and outputs a fasta file on stdin.
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('contaminantFile',
                    help='Filepath to contaminant list file')

opts = parser.parse_args()

f = open(opts.contaminantFile,'r')
for line in f:
    line = line.strip()
    if(line.startswith('#')):
        continue
    sline = line.split('\t')
    if(len(sline)<2):
        continue
    header = sline[0]
    seq = sline[-1]
    print(">%s\n%s" % (header,seq))

f.close()
