from seqtools.snpeff import effNames,snpEffEffects

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
    header += ['FILTER', 'CHROM', 'POS', 'REF', 'ALT', 'ID'] + ['info.' + x for x in infos]
    if(snpeff):
        header += effNames
    formatList = []
    for format in formats:
        for sample in reader.samples:
            formatList.append(sample + "." + format)
    header += formatList 

    outfile.write('\t'.join(header) + "\n")
        


    def flatten(x):
        if type(x) == type([]):
            # can probably change to tab-separated if headers match up right!
            x = ','.join(map(str, x))
        return x

    for record in reader:
        info_row = [flatten(record.INFO.get(x, None)) for x in infos]
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
                row.append(flatten(getattr(sample.data, x, '')))
        newrow = []
        for r in row:
            if(r is None):
                newrow.append('')
            else:
                newrow.append(str(r))
        outfile.write('\t'.join(newrow) + "\n")
