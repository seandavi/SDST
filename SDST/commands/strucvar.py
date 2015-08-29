from SDST.main import subparsers
import SDST.strucvar
import SDST.strucvar.crest

def crestToBed(opts):
    with open(opts.predSVfile,'r') as f:
        for line in f:
            if(opts.samplename is not None):
                print(SDST.strucvar.crest.crestLineToBedLines(line,opts.samplename))
            else:
                print(SDST.strucvar.crest.crestLineToBedLines(line))
                


strucvar_parser = subparsers.add_parser('strucvar')
strucvar_subparsers = strucvar_parser.add_subparsers(help="Structural variation subcommands")
crestToBed_parser = strucvar_subparsers.add_parser('crestToBed',
                                                   help="convert CREST output to two BED records per line")
crestToBed_parser.add_argument('-s','--samplename',default=None,type=str,
                               help="sample information to add to the BED 'name' column")
crestToBed_parser.add_argument('predSVfile',help="The name of the CREST predSV.txt file")

crestToBed_parser.set_defaults(func=crestToBed)

