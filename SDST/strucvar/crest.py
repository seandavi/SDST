def crestLineToBedLines(crestline,extrastring=None):
    """Takes a line from a CREST file and turns it into two BED lines

    :param crestline: a single string representing the CREST output
    :param extrastring: a single string to concatenate to the bed output, useful for including sample information, etc.
    :returns: a string containing the two bed lines"""
    sline = crestline.strip().split("\t")
    outstringparts=[sline[8],sline[0],sline[1],sline[4],sline[5]]
    if(extrastring is not None):
        outstringparts+=[str(extrastring)]
    namefield = "_".join(outstringparts)
    leftbreak = "\t".join([sline[0],sline[1],sline[1],namefield])
    rightbreak = "\t".join([sline[4],sline[5],sline[5],namefield])
    return(leftbreak+"\n"+rightbreak)
