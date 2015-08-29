#!/usr/bin/env python
"""
parse VCF output from VarScan
and fix the ALT column
to adhere with VCF specifications
"""
import re

def fixLine(line):
    """Fix a varscan VCF line

    Prints the output to stdout.  Fixes the ALT column and also fixes the FREQ field to be a floating point value, easier for filtering.

    :param line: a pre-split and stripped varscan line
    """
    line = line.strip()
    
    if(line.startswith("##")):
        line=line.replace('FREQ,Number=1,Type=String',
                          'FREQ,Number=1,Type=Float')
        return line
    
    if(line.startswith("#CHROM")):
        return line

    line = line.split('\t')
    try:
        REF,ALT = line[3:5]
    except ValueError:
        return "\t".join(line) + "\n"

    Ifreq = line[8].split(":").index("FREQ")
    ndat = line[9].split(":")
    tdat = line[10].split(":")
    ndat[Ifreq] = str(float(ndat[Ifreq].rstrip("%"))/100)
    tdat[Ifreq] = str(float(tdat[Ifreq].rstrip("%"))/100)
    line[9]=":".join(ndat)
    line[10]=":".join(tdat)
    
    
    if "+" in ALT or "-" in ALT:
        if "/" not in ALT:
            if ALT[0] == "+":
                R = REF
                A = REF + re.sub(r'^[+-][\d]?','',ALT)
            elif ALT[0] == "-":
                R = REF + re.sub(r'^[+-][\d]?','',ALT)
                A = REF
        else:
            Ins = [p[1:] for p in ALT.split("/") if p[0]=="+"]
            Del = [p[1:] for p in ALT.split("/") if p[0]=="-"]

            if len(Del):
                REF += sorted(Del,key= lambda x: len(x))[-1]

            A = ",".join([REF[::-1].replace(p[::-1], "", 1)[::-1] for p in Del] + [REF+p for p in Ins])
            R = REF

        REF = R
        ALT = A
    else:
        ALT = ALT.replace('/',',')

    line[3] = REF
    line[4] = ALT
    return "\t".join(line)

def fixVarscanVcfFile(iterable):
    """Takes an interator over a varscan VCF file and returns an iterator over fixed VCF lines, including header.
    
    :param iterable: any iterable of the VCF lines
    :returns: An iterator over fixed VCF lines

    Usage is like so:

    >>> from SDST.varscan import fixVarscanVcfFile
    >>> varscan = fixVarscanVcfFile(open('filename.vcf','r'))
    >>> for line in varscan:
        print line
        
    """
    for line in iterable:
        yield fixLine(line)

def main():
    import argparse,sys
    parser = argparse.ArgumentParser('Parse VCF output from Varscan to output valid VCF.  Output is to stdout.')

    parser.add_argument('-v','--varscan',
                        help="varscan vcf output file name")

    opts = parser.parse_args()
    if(not opts.varscan):
        varscan = fixVarscanVcfFile(sys.stdin)
    else:
        varscan = fixVarscanVcfFile(open(opts.varscan))
    for line in varscan:
        print(line)

if __name__ == '__main__':
    main()

