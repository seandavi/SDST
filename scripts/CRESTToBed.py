#!/usr/bin/env python
import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument('predSVfile')
args = parser.parse_args()

fname = os.path.basename(args.predSVfile)

with open(args.predSVfile,'r') as f:
    for line in f:
        sline = line.strip().split("\t")
        name = "_".join([sline[8],sline[0],sline[1],sline[4],sline[5],fname])
        leftbreak = "\t".join([sline[0],sline[1],sline[1],name])
        rightbreak = "\t".join([sline[4],sline[5],sline[5],name])
        print(leftbreak)
        print(rightbreak)
