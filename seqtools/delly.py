class DellyParser(object):
    """Return an iterator of DellyRecords
    """
    
    def __init__(self, filename):
        """

        Arguments:
        - `filename`: file name of the delly output file
        """
        
        

class DellyReadRecord(object):
    """Represents a read record for a delly record
    """
    
    def __init__(self, line):
        """Initialize with a simple string representing the read record
        
        Arguments:
        - `line`:
        """
        self._line = line
        sline = line.strip().split()
        self.rname       = sline[0]
        self.filter      = sline[1]
        self.left_chrom  = sline[2]
        self.left_pos    = int(sline[3])
        self.mapq        = int(sline[4])
        self.right_chrom = sline[5]
        self.right_pos   = int(sline[6])
        self.other       = int(sline[7])
        self.library     = sline[8]
        
        

class DellyRecord(object):
    """Encapsulate an entire delly record including the reads
    """
    
    def __init__(self):
        """
        """
        self._reads = []

    def addRead(line):
        """
        Add a read record to the DellyRecord
        `line`: the line as a string to add
        """
        self._reads.append(DellyReadRecord(line))
