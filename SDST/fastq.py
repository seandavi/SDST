from SDST.utils import fileOpen

class FastqRecord(object):
    """Very simple fastq class containing header, sequence, line3, and quality as strings"""
    def __init__(self,header,sequence,line3,quality):
        """Initialize a FastqRecord

        :param header: header
        :param sequence: sequence
        :param line3: second header line from fastq record
        :param quality: string quality"""
        self.header=header
        self.sequence=sequence
        self.line3=line3
        self.quality=quality

    def __str__(self):
        """String represendation, suitable for output to a fastq file"""
        return "{}\n{}\n{}\n{}".format(self.header,self.sequence,self.line3,self.quality)
    

class Fastq(object):
    def __init__(self,fname):
        self.name = fname
        self.fh = fileOpen(fname)
        
    def _getNextRecord(self):
        x = []
        for i in range(0,4):
            nextLine=self.fh.readline().strip()
            x.append(nextLine)
            if(nextLine==''):
                return None
        return(FastqRecord(x[0],x[1],x[2],x[3]))

    def __iter__(self):
        return self

    def __next__(self):
        tmp = self._getNextRecord()
        if(tmp is None):
            raise StopIteration
        else:
            return(tmp)

    def name(self):
        return self.fname

    def close(self):
        self.fh.close()
