import subprocess
import tempfile
import gzip

try:
    transtab = str.maketrans('ACGTNacgtn','TGCANtgcan')
except:
    import string
    transtab = string.maketrans('ACGTNacgtn','TGCANtgcan')


revcompdict = {}
def revcomp(sequence):
    """Reverse complement a string

    :param sequence: The DNA string (all caps)

    This function includes a caching mechanism, so watch memory usage!
    """
    if(sequence in revcompdict):
        return revcompdict[sequence]
    tmp=sequence[::-1].translate(transtab)
    revcompdict[sequence]=tmp
    return(tmp)


def fileOpen(fname,mode='rt',encoding='latin-1'):
    """Open a file, including gzip files

    :param fname: The filename to open.  Gzip files are distinguished by ending in '.gz'
    :param mode: File mode for opening [default 'r']
    :returns: An open file handle

    Note: the gzip module in python is REALLY slow, so this function uses
    subprocess and command-line gzip instead.
    """
    # gzip in python is REALLY slow, so use pipes instead.
    if(fname.endswith('.gz')):
        return gzip.open(fname,mode=mode,encoding=encoding)
    else:
        return open(fname,mode)


def sortVcfBySequence(vcf,seqnames,seqmap=None):
    """Sort a tabixed VCF file by sequence names

    >>> import vcf
    >>> v = vcf.Reader('vcf.gz')
    >>> from SDST.utils import sortVcfBySequence
    >>> w = vcf.Writer(open('sorted.vcf','w'),v)
    >>> faifile = open('ucsc.hg19.fasta.fai')
    >>> seqnames = [x.split('\t')[0] for x in faifile]
    >>> for rec in sortVcfBySequence(v,seqnames):
        w.write_record(rec)
    >>> w.close()
    >>> v.close()"""

    if(seqmap):
        revseqmap = dict([(v,k) for (k,v) in list(seqmap.items())])

    print(seqmap)
    
    for s in seqnames:
        # seqmap maps seqnames to tabix index names
        if(seqmap is not None):
            region = vcf.fetch(seqmap[s],start=0)
        else:
            try:
                region = vcf.fetch(s,start=0)
                print(region)
            except ValueError:
                print("here",s)
                continue
        print(region)
        print(s)
        for rec in region:
            if(seqmap is not None):
                rec.CHROM=revseqmap[rec.CHROM]
            yield rec
    

