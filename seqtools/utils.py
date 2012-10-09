from string import maketrans
import subprocess

transtab = maketrans('ACGTNacgtn','TGCANtgcan')

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


def fileOpen(fname,mode='r'):
    """Open a file, including gzip files

    :param fname: The filename to open.  Gzip files are distinguished by ending in '.gz'
    :param mode: File mode for opening [default 'r']
    :returns: An open file handle

    Note: the gzip module in python is REALLY slow, so this function uses
    subprocess and command-line gzip instead.
    """
    # gzip in python is REALLY slow, so use pipes instead.
    if(fname.endswith('.gz')):
        if(mode.startswith('r')):
            return subprocess.Popen(['gunzip -c %s' % fname],stdout=subprocess.PIPE,shell=True).stdout
        if(mode.startswith('w')):
            return subprocess.Popen(['gzip > %s' % fname],stdin=subprocess.PIPE,shell=True).stdin
    else:
        return open(fname,'r')

