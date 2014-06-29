#!/usr/bin/env python
import vcf
import argparse
import sys

parser = argparse.ArgumentParser(description='Renames the samples in a VCF file taken from stdin and writes out to stdout.')

parser.add_argument('-s','--sampleName',action='append',
                    help='Specify names, one name per instance of -n ("-n sample1 -n sample2 -n sample3")')
opts = parser.parse_args()

reader = vcf.Reader(sys.stdin)
if(len(reader.samples)!=len(opts.sampleName)):
    sys.stderr.write('The number of names specified needs to match the number of samples in the VCF files')
    exit(-1)
reader.samples=opts.sampleName
writer = vcf.Writer(sys.stdout,reader)
for record in reader:
    # sys.stderr.write(str(record) + "\n")
    writer.write_record(record)

writer.close()
