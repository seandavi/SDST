import vcf.parser

def addStrelkaHeaders(reader):
    """Add strelka headers to a VCF reader
    
    :param reader: a VCR reader object from pyVCF
    
    The headers 'TUMREF','TUMALT','NORMREF', and 'NORMALT' are added.

    After the reader is updated with new information, the records 
    can be manipulated to add new information and the writer will
    know what to do with it."""

    reader.infos['TUMREF'] = vcf.parser._Info(id='TUMREF',num=1,type='Integer',desc="The depth of the REF allele in the tumor")
    reader.infos['TUMALT'] = vcf.parser._Info(id='TUMALT',num=1,type='Integer',desc="The depth of the FIRST ALT allele in the tumor")
    reader.infos['NORMREF'] = vcf.parser._Info(id='NORMREF',num=1,type='Integer',desc="The depth of the REF allele in the normal")
    reader.infos['NORMALT'] = vcf.parser._Info(id='NORMALT',num=1,type='Integer',desc="The depth of the FIRST ALT allele in the normal")
    reader.infos['TUMVAF'] = vcf.parser._Info(id='TUMVAF',num=1,type='Float',desc="The Variant Allele Frequency in the tumor")
    reader.infos['TUMVARFRACTION'] = vcf.parser._Info(id='TUMVARFRACTION',num=1,type='Float',desc="The fraction of variant reads in the tumor versus the total (so 1.0 means all variant reads in the tumor)")
    
def modifyStrelkaRow(record,fixIndels=True):
    """Add info for strelka processing to vcf record
    
    :param record: a pyVCF record object
    """
    if(record.is_snp or record.ALT[0] is None):
        ref = record.REF
        alt = record.ALT[0]
        record.INFO['NORMREF']=getattr(record.samples[0].data,ref+'U')[0]
        record.INFO['TUMREF']=getattr(record.samples[1].data,ref+'U')[0]
        # strelka sometimes reports a non-passing variant as no "ALT" allele (no change)
        if(alt is None): 
            record.INFO['NORMALT']=0
            record.INFO['TUMALT']=0
            record.INFO['TUMVAF']=0
            record.INFO['TUMVARFRACTION']=0
        else:
            alt = str(alt)
            record.INFO['NORMALT']=getattr(record.samples[0].data,alt+'U')[0]
            record.INFO['TUMALT']=getattr(record.samples[1].data,alt+'U')[0]
            try:
                record.INFO['TUMVAF']=float(record.INFO['TUMALT'])/(record.INFO['TUMALT']+record.INFO['TUMREF'])
            except ZeroDivisionError:
                record.INFO['TUMVAF']=0
            try:
                record.INFO['TUMVARFRACTION']=float(record.INFO['TUMALT'])/(record.INFO['TUMALT']+record.INFO['NORMALT'])
            except ZeroDivisionError:
                record.INFO['TUMVARFRACTION']=0
        return(record)
    else:
        record.INFO['NORMREF']=getattr(record.samples[0].data,'TAR')[0]
        record.INFO['NORMALT']=getattr(record.samples[0].data,'TIR')[0]
        record.INFO['TUMREF']=getattr(record.samples[1].data,'TAR')[0]
        record.INFO['TUMALT']=getattr(record.samples[1].data,'TIR')[0]
        try:
            record.INFO['TUMVAF']=float(record.INFO['TUMALT'])/(record.INFO['TUMALT']+record.INFO['TUMREF'])
        except ZeroDivisionError:
            record.INFO['TUMVAF']=0
        try:
            record.INFO['TUMVARFRACTION']=float(record.INFO['TUMALT'])/(record.INFO['TUMALT']+record.INFO['NORMALT'])
        except ZeroDivisionError:
            record.INFO['TUMVARFRACTION']=0
        if(fixIndels):
            record.REF=record.REF.replace('.','')
            for i in range(len(record.ALT)):
                if(not isinstance(record.ALT[i],vcf.model._Substitution)):
                    return(None)
        return(record)
        

def processVcf(reader,outfile):
    addStrelkaHeaders(reader)
    writer = vcf.Writer(outfile,reader)
    for row in reader:
        newrow = modifyStrelkaRow(row)
        if(newrow is not None):
            writer.write_record(newrow)
