from seqtools.snpeff import effNames,snpEffEffects
import itertools
import vcf

def vcfMelt(reader,outfile,samplename=None):
    """Melt a VCF file into a tab-delimited text file

    :param reader: a vcf.Reader object (from the pyvcf package)
    :param outfile: a stream to which to write the resulting melted VCF file
    """
    formats = list(reader.formats.keys())
    infos = list(reader.infos.keys())
    if('EFF' in reader.infos):
        snpeff = True
        del infos[infos.index('EFF')]
    # TODO: split out formats per sample
    # TODO: split out filters per column
    # TODO: split up info fields and format array fields into separate columns
    header = []
    if(samplename is not None):
        header += ['SampleName']
    header += ['FILTER', 'CHROM', 'POS', 'REF', 'ALT', 'ID'] 
    for x in infos:
        infonum=reader.infos.get(x)[1]
        if(infonum is None): infonum=0
        if(infonum>1):
            for y in range(infonum):
                header.append(x + "." + str(y))
        else:
            header.append(x)
    if(snpeff):
        header += effNames
    formatList = []
    for format in formats:
        for sample in reader.samples:
            infonum=reader.formats.get(format)[1]
            if(infonum is None): infonum=0
            if(infonum>1):
                for y in range(infonum):
                    formatList.append(sample + "." + format + "." + str(y))
            else:
                formatList.append(sample + "." + format)
    header += formatList 

    outfile.write('\t'.join(header) + "\n")


    def flatten(x):
        if type(x) == type([]):
            # can probably change to tab-separated if headers match up right!
            return(list(map(str, x)))
        else:
            return([str(x)])

    for record in reader:
        info_row=[]
        for x in infos:
            infonum = reader.infos.get(x)[1]
            if(infonum is None):
                infonum = 0
            val = record.INFO.get(x,None)
            if((type(val)!=type([])) and (infonum>1)):
                val = list(itertools.repeat('',infonum))
            if((infonum<=1) and (type(val)==type([]))):
                val = ','.join(str(v) for v in val)
            info_row += flatten(val) 
        if(snpeff):
            try:
                maxeffect = snpEffEffects(record.INFO['EFF']).highest
                info_row += list(maxeffect.values())
            except KeyError:
                # return a bunch of NAs
                info_row += ['NA']*len(effNames)
        fixed = [record.CHROM, record.POS, record.REF, ','.join([str(alt) for alt in record.ALT]), record.ID] 

        row = []
        if(samplename is not None):
            row += [samplename]
        if(record.FILTER):
            row += [",".join(record.FILTER)]
        else:
            row += ['.']
        row += fixed
        row += info_row

        for x in formats:
            for sample in record.samples:
                row+=flatten(getattr(sample.data, x, ''))
        newrow=[]
        for r in row:
            if(r is None):
                newrow.append('')
            else:
                newrow.append(str(r))
        outfile.write('\t'.join(newrow) + "\n")


def basesAtPos(samfile, pos, chromname, minbasequal, minmapqual):
    'Return a string of the bases at that position.'
    position = 0
    coverage = 0
    bases = ""
    for pileupcolumn in samfile.pileup(reference=chromname, start=pos-1, end=pos):
        if ((pileupcolumn.pos+1)==pos):
            position = int(pileupcolumn.pos+1)
            coverage = int(pileupcolumn.n)
            for pileupread in pileupcolumn.pileups:
                if(pos==73433494):
                    print(pos,bases)
                if (pileupread.indel == 0 and pileupread.is_del == 0 and \
                (pileupread.alignment.qual[pileupread.qpos]) >= minbasequal and \
                float(pileupread.alignment.mapq) >= minmapqual):
                    bases += chr(pileupread.alignment.seq[pileupread.qpos])
    return position, coverage, bases


def countBases(reader,outfile,bamfile):
    """Count the ref and alt bases in a bam file at each variant location

    :params reader: a VCF Reader object
    :params outfile: a stream to which to write the output
    :params bamfile: a pysam Samfile, sorted and indexed

    Adds the INFO fields:

    RNAC_REF, RNAC_ALT, RNAC_MAF"""

    reader.infos['RNAC_REF'] = vcf.parser._Info(id='RNAC_REF',num=1,type='Integer',desc='The count of REF alleles in the bamfile')    
    reader.infos['RNAC_ALT'] = vcf.parser._Info(id='RNAC_ALT',num=1,type='Integer',desc='The count of ALT alleles in the bamfile')    
    reader.infos['RNAC_MAF'] = vcf.parser._Info(id='RNAC_MAF',num=1,type='Float',desc='The fraction of ALT allele in the bamfile')

    writer = vcf.Writer(outfile,reader)

    for row in reader:
        ref = row.REF
        alt = str(row.ALT[0])
        bases = basesAtPos(bamfile,row.POS,row.CHROM,0,0)[2]
        refcount,altcount = [len([x for x in bases if x==ref]),
                             len([x for x in bases if x==alt])]
        row.INFO['RNAC_REF']=refcount
        row.INFO['RNAC_ALT']=altcount
        row.INFO['RNAC_MAF']=0
        if(refcount+altcount>0):
            row.INFO['RNAC_MAF']=float(altcount)/(refcount+altcount)

        writer.write_record(row)


