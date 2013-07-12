import sys
import argparse
from seqtools.main import subparsers
import seqtools.vcf
import seqtools.strelka
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
    if(opts.includeHeader):
        writer = vcf.Writer(outfile,v)
    seqtools.vcf.vcfMelt(v,outfile,opts.samplename)
    outfile.close()

def strelkaProcess(opts):
    if(opts.vcf):
        v = vcf.Reader(filename=opts.vcf)
    else:
        v = vcf.Reader(sys.stdin)
    if(opts.outfile):
        outfile = open(opts.outfile,'w')
    else:
        outfile = sys.stdout
    seqtools.strelka.processVcf(v,outfile)
    outfile.close()


varscan_parser = subparsers.add_parser('vcf')
varscan_subparsers = varscan_parser.add_subparsers(help="vcf subcommands")
meltVcf_parser = varscan_subparsers.add_parser('melt',
                                              help="Melt a VCF file into a tab-delimited text file with special treatment of snpEff annotation")
meltVcf_parser.add_argument('-f','--vcf',
                            help='Filename of VCF file [default=stdin]')
meltVcf_parser.add_argument('-o','--outfile',
                            help='Filename of VCF file [default=stdout]')
meltVcf_parser.add_argument('-s','--samplename',default=None,
                            help='Sample name to include in first column of output [default=None]')
meltVcf_parser.add_argument('-i','--includeHeader',action='store_true',
                            help='Include VCF header in the output; useful for including definitions of columns')

strelka_parser = varscan_subparsers.add_parser('strelka',
                                              help="process a vcf strelka-produced vcf file to add tumor/normal read count info fields (TUMREF,NORMREF,TUMALT,NORMALT,TUMVAF)")
strelka_parser.add_argument('-f','--vcf',
                            help='Filename of VCF file [default=stdin]')
strelka_parser.add_argument('-o','--outfile',
                            help='Filename of VCF file [default=stdout]')

strelka_parser.set_defaults(func=strelkaProcess)
