__author__ = 'sdavis2'

import vcf
import argparse

parser = argparse.ArgumentParser(description="""Convert all structural variant lines into bed lines""")

parser.add_argument('-s','--minsomaticscore',default=15,type=int)
parser.add_argument("vcffile")

opts = parser.parse_args()

def strBND(line):
    """
    :param line: a VCF Record representing a Breakend
    :return: A simple bed format line (stringified, already)
    """
    return("\t".join([str(x) for x in (line.CHROM,line.POS,line.POS+1,line.ALT[0])]))
with open(opts.vcffile) as f:
    reader = vcf.Reader(f)
    for line in reader:
        if((line.is_sv)):
            if(line.INFO['SOMATICSCORE']<opts.minsomaticscore):
                continue
            if(line.INFO['SVTYPE']=="BND"):
                print(strBND(line))
            else:
                print("\t".join([str(x) for x in (line.CHROM,line.POS,line.INFO['END'],line.ALT[0])]))
