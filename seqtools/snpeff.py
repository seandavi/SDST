import re
import collections 

impactOrder = ('HIGH','MODERATE','LOW','MODIFIER')

effNames = ['Effect',
            'Impact',
            'Functional',
            'CodonChange',
            'AAChange',
            'AALength',
            'Gene',
            'BioType',
            'Coding',
            'Transcript',
            'Exon',
            'Warnings']

class snpEffRecord(collections.OrderedDict):
    def __init__(self,*effs):
        super(snpEffRecord,self).__init__(zip(effNames,effs))



effRe = re.compile('(\w+)\((.*)\)')

class snpEffEffects(object):
    """Operate on a snpEff VCF string
    """
    def __init__(self,eff):
        tmp = eff.split(',')
        self.effs = []
        for effgroup in tmp:
            (maineff,mods) = effRe.match(effgroup).groups()
            efflist = [maineff] + mods.split('|')
            if(len(efflist)==11):
                efflist.append('')
            self.effs.append(snpEffRecord(*efflist))

    @property
    def highest(self,impactOrder=impactOrder):
        for impact in impactOrder:
            for eff in self.effs:
                if(eff['Impact'] == impact):
                    return eff
        return None

    
            
