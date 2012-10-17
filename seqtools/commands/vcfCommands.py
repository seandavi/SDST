import sys
import argparse
from seqtools.main import subparsers
import seqtools.vcf
import vcf

def meltVcf(opts):
    if(opts.vcf):
        v = vcf.Reader(filename=opts.vcf)
    else:
        v = vcf.Reader(sys.stdin)
    if(opts.outfile):
        outfile = open(opts.outfile,'w')
    else:
        outfile = sys.stdout
    print v,outfile
    seqtools.vcf.vcfMelt(v,outfile)
    outfile.close()

varscan_parser = subparsers.add_parser('vcf')
varscan_subparsers = varscan_parser.add_subparsers(help="vcf subcommands")
meltVcf_parser = varscan_subparsers.add_parser('melt',
                                              help="Melt a VCF file into a tab-delimited text file with special treatment of snpEff annotation")
meltVcf_parser.add_argument('-f','--vcf',
                            help='Filename of VCF file [default=stdin]')
meltVcf_parser.add_argument('-o','--outfile',
                            help='Filename of VCF file [default=stdout]')
meltVcf_parser.set_defaults(func=meltVcf)
