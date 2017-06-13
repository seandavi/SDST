import sys
import argparse
from SDST.main import subparsers
import SDST.seqvcf
import SDST.strelka
import SDST.mutect
import vcf
import pysam

def meltVcf(opts):
    if(opts.vcf):
        f = open(opts.vcf,mode='r')
        v = vcf.Reader(f)
    else:
        v = vcf.Reader(sys.stdin)
    if(opts.outfile):
        outfile = open(opts.outfile,'w')
    else:
        outfile = sys.stdout
    if(opts.includeHeader):
        writer = vcf.Writer(outfile,v)
    SDST.seqvcf.vcfMelt(v,outfile,opts.samplename,opts.includeGenotypes)
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
    SDST.strelka.processVcf(v,outfile)
    outfile.close()

def mutectProcess(opts):
    if(opts.vcf):
        v = vcf.Reader(filename=opts.vcf)
    else:
        v = vcf.Reader(sys.stdin)
    if(opts.outfile):
        outfile = open(opts.outfile,'w')
    else:
        outfile = sys.stdout
    SDST.mutect.processVcf(v,outfile)
    outfile.close()

def RNAcounts(opts):
    if(opts.vcf):
        f = open(opts.vcf,mode='r')
        v = vcf.Reader(f)
    else:
        v = vcf.Reader(sys.stdin)
    if(opts.outfile):
        outfile = open(opts.outfile,'w')
    else:
        outfile = sys.stdout
    bamfile = pysam.Samfile(opts.bamfile)
    SDST.seqvcf.countBases(v,outfile,bamfile)
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
meltVcf_parser.add_argument('-g','--includeGenotypes',action='store_true',
                            help='Include genotypes for each sample (XXX.GT in column header) if available')
meltVcf_parser.set_defaults(func=meltVcf)


strelka_parser = varscan_subparsers.add_parser('strelka',
                                              help="process a vcf strelka-produced vcf file to add tumor/normal read count info fields (TUMREF,NORMREF,TUMALT,NORMALT,TUMVAF)")
strelka_parser.add_argument('-f','--vcf',
                            help='Filename of VCF file [default=stdin]')
strelka_parser.add_argument('-o','--outfile',
                            help='Filename of VCF file [default=stdout]')

strelka_parser.set_defaults(func=strelkaProcess)


mutect_parser = varscan_subparsers.add_parser('mutect',
                                              help="process a vcf mutect-produced vcf file to add tumor/normal read count info fields (TUMREF,NORMREF,TUMALT,NORMALT,TUMVAF)")
mutect_parser.add_argument('-f','--vcf',
                            help='Filename of VCF file [default=stdin]')
mutect_parser.add_argument('-o','--outfile',
                            help='Filename of VCF file [default=stdout]')

mutect_parser.set_defaults(func=mutectProcess)




rnacount_parser = varscan_subparsers.add_parser('rnacount',
                                                help="Count the number of REF and ALT alleles in a bamfile at each variant location in a VCF file.  The results are added to the VCF file in three new INFO tags, RNAC_REF, RNAC_ALT, and RNAC_MAF")

rnacount_parser.add_argument('-f','--vcf',
                            help='Filename of VCF file [default=stdin]')
rnacount_parser.add_argument('-o','--outfile',
                            help='Filename of VCF file [default=stdout]')
rnacount_parser.add_argument('bamfile',
                             help='The name of the RNA-seq bamfile from which to collect counts')

rnacount_parser.set_defaults(func=RNAcounts)










