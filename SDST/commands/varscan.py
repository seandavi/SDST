import sys
import argparse
from SDST.main import subparsers
from SDST.varscan import fixVarscanVcfFile


def fixVcf(opts):
    if(not opts.varscan):
        varscan = fixVarscanVcfFile(sys.stdin)
    else:
        varscan = fixVarscanVcfFile(opts.varscan)
    for line in varscan:
        print(line)

varscan_parser = subparsers.add_parser('varscan')
varscan_subparsers = varscan_parser.add_subparsers(help="VarScan subcommands")
fixVcf_parser = varscan_subparsers.add_parser('fixVcf',
                                              help="Fix the FREQ field and the ALT field in varscan VCF output")
fixVcf_parser.add_argument('-f','--varscan',default=sys.stdin,
                           type=argparse.FileType('r'),
                           help='Filename of varscan VCF file [default=stdin]')
fixVcf_parser.set_defaults(func=fixVcf)
