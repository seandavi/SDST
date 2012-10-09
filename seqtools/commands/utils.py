from seqtools.main import subparsers
from seqtools.demultiplexer import demultiplex

def demultiplex_setup(opts):
    if(opts.readFile2):
        demultiplex(readfile = opts.readFile1,
                    readfile2 = opts.readFile2,
                    indexfile = opts.indexFile,
                    indexes = opts.index)
    else:
        demultiplex(readfile = opts.readFile1,
                    indexfile = opts.indexFile,
                    indexes = opts.index)

utils_parser = subparsers.add_parser('utils')
utils_subparsers = utils_parser.add_subparsers(help="General utilities")
demultiplex_parser = utils_subparsers.add_parser('demultiplex',
                                                 help="Demultiplex a set of fastq files")
demultiplex_parser.add_argument('-1','--readFile1',
                    help="read1 filename")
demultiplex_parser.add_argument('-2','--readFile2',
                    help="read2 filename")
demultiplex_parser.add_argument('-i','--indexFile',
                    help="index Filename")
demultiplex_parser.add_argument('-x','--index',type=str,action='append',
                    help="The indexes, one per index")

demultiplex_parser.set_defaults(func=demultiplex_setup)
