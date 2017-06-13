import vcf
import fisher
import math

def addMutectHeaders(reader):
    """Add mutect headers to a VCF reader
    
    :param reader: a VCF reader object from pyVCF
    
    The headers 'TUMREF','TUMALT','NORMREF', and 'NORMALT' are added.

    After the reader is updated with new information, the records 
    can be manipulated to add new information and the writer will
    know what to do with it."""

    reader.infos['TUMREF'] = vcf.parser._Info(id='TUMREF',num=1,type='Integer',desc="The depth of the REF allele in the tumor",source=None,version=None)
    reader.infos['TUMALT'] = vcf.parser._Info(id='TUMALT',num=1,type='Integer',desc="The depth of the FIRST ALT allele in the tumor",source=None,version=None)
    reader.infos['NORMREF'] = vcf.parser._Info(id='NORMREF',num=1,type='Integer',desc="The depth of the REF allele in the normal",source=None,version=None)
    reader.infos['NORMALT'] = vcf.parser._Info(id='NORMALT',num=1,type='Integer',desc="The depth of the FIRST ALT allele in the normal",source=None,version=None)
    reader.infos['TUMVAF'] = vcf.parser._Info(id='TUMVAF',num=1,type='Float',desc="The Variant Allele Frequency in the tumor",source=None,version=None)
    reader.infos['TUMVARFRACTION'] = vcf.parser._Info(id='TUMVARFRACTION',num=1,type='Float',desc="The fraction of variant reads in the tumor versus the total (so 1.0 means all variant reads in the tumor)",source=None,version=None)
    reader.infos['LOG_FISHER'] = vcf.parser._Info(id='LOG_FISHER',num=1,type='Float',desc="-Log10(Fisher Exact Test p-value)",source=None,version=None)
    
def modifyMutectRow(record,fixIndels=True):
    """Add info for mutect processing to vcf record
    
    :param record: a pyVCF record object
    """
    record.INFO['NORMREF']=record.genotype('NORMAL')['AD'][0]
    record.INFO['TUMREF']=record.genotype('TUMOR')['AD'][0]
    if(record.ALT is None): 
        record.INFO['NORMALT']=0
        record.INFO['TUMALT']=0
        record.INFO['TUMVAF']=0
        record.INFO['TUMVARFRACTION']=0
    else:
        record.INFO['NORMALT']=record.genotype('NORMAL')['AD'][1]
        record.INFO['TUMALT']=record.genotype('TUMOR')['AD'][1]
        try:
            record.INFO['TUMVAF']=float(record.INFO['TUMALT'])/(record.INFO['TUMALT']+record.INFO['TUMREF'])
        except ZeroDivisionError:
            record.INFO['TUMVAF']=0
        try:
            record.INFO['TUMVARFRACTION']=float(record.INFO['TUMALT'])/(record.INFO['TUMALT']+record.INFO['NORMALT'])
        except ZeroDivisionError:
            record.INFO['TUMVARFRACTION']=0
        record.INFO['LOG_FISHER']=-math.log10(fisher.pvalue(record.INFO['TUMREF'],record.INFO['TUMALT'],record.INFO['NORMREF'],record.INFO['NORMALT']).two_tail)
    return(record)

def processVcf(reader,outfile):
    addMutectHeaders(reader)
    writer = vcf.Writer(outfile,reader)
    for row in reader:
        newrow = modifyMutectRow(row)
        if(newrow is not None):
            writer.write_record(newrow)
