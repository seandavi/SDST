from SDST.snpeff import effNames,snpEffEffects
import itertools
import vcf
import math
import fisher

def vcfMelt(reader,outfile,samplename=None,includeGenotypes=False):
    """Melt a VCF file into a tab-delimited text file

    :param reader: a vcf.Reader object (from the pyvcf package)
    :param outfile: a stream to which to write the resulting melted VCF file
    :param samplename: include a samplename as the first column for the case for the case of multiple text files being concatenated in "long" format -- default "None" for not including
    :param includeGenotypes: Include genotypes in output (one per sample in the XXX.GT columns) -- default False
    """
    formats = list(reader.formats.keys())
    infos = list(reader.infos.keys())
    snpeff = False
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
    if(includeGenotypes):
        for sample in reader.samples:
            formatList.append(sample + ".GT")
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
                if(maxeffect is None):
                    # effect is as a modifier and chromosome is not found
                    info_row += ['NA']*len(effNames)
                else:
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

        if(includeGenotypes):
            for sample in reader.samples:
                calldata=record.genotype(sample)
                if(calldata.called):
                    row.append(str(record.genotype(sample).data.GT))
                else:
                    row.append('./.')
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
                (pileupread.alignment.query_qualities[pileupread.query_position]) >= minbasequal and \
                float(pileupread.alignment.mapping_quality) >= minmapqual):
                    bases += pileupread.alignment.query_sequence[pileupread.query_position]
    return position, coverage, bases


def countBases(reader,outfile,bamfile,prefix='RNAC',vaf='MAF'):
    """Count the ref and alt bases in a bam file at each variant location

    :params reader: a VCF Reader object
    :params outfile: a stream to which to write the output
    :params bamfile: a pysam Samfile, sorted and indexed
    :params prefix: The prefix string in the INFO fields that will be added to the VCF file
    :params vaf: The suffix string added for the variant allele frequency name in the VCF INFO field

    Adds the INFO fields:

    {prefix}_REF, {prefix}_ALT, {prefix}_{vaf}"""

    REFCOUNT_NAME = prefix + "_REF"
    ALTCOUNT_NAME = prefix + "_ALT"
    VAF_NAME      = prefix + "_" + vaf

    reader.infos[REFCOUNT_NAME] = vcf.parser._Info(id=REFCOUNT_NAME,num=1,type='Integer',desc='The count of REF alleles in the bamfile',source='',version='')
    reader.infos[ALTCOUNT_NAME] = vcf.parser._Info(id=ALTCOUNT_NAME,num=1,type='Integer',desc='The count of ALT alleles in the bamfile',source='',version='')
    reader.infos[VAF_NAME     ] = vcf.parser._Info(id=VAF_NAME,num=1,type='Float',desc='The fraction of ALT allele in the bamfile',source='',version='')

    writer = vcf.Writer(outfile,reader)

    for row in reader:
        ref = row.REF
        alt = str(row.ALT[0])
        bases = basesAtPos(bamfile,row.POS,row.CHROM,0,0)[2]
        refcount,altcount = [len([x for x in bases if x==ref]),
                             len([x for x in bases if x==alt])]
        row.INFO[REFCOUNT_NAME]=refcount
        row.INFO[ALTCOUNT_NAME]=altcount
        row.INFO[VAF_NAME]=0
        if(refcount+altcount>0):
            row.INFO[VAF_NAME]=float(altcount)/(refcount+altcount)

        writer.write_record(row)

def addFisher(reader,outfile):
    writer = vcf.Writer(outfile,reader)
    reader.infos['LOG_FISHER'] = vcf.parser._Info(id='LOG_FISHER',num=1,type='Float',desc='Log10 of Fisher Exact Test p-value',source='',version='')
    
    for row in reader:
        tumorcounts = row.genotype('TUMOR')['AD']
        normalcounts = row.genotype('NORMAL')['AD']
        row.INFO['LOG_FISHER'] = -math.log10(fisher.pvalue(tumorcounts[0],tumorcounts[1],normalcounts[0],normalcounts[1]).two_tail)
        writer.write_record(row)
