#!/usr/bin/env python
CLSPRIORITY = [
    'frameshift_element', 'thermoregulator', 'IRES', 'riboswitch', 'ribozyme',
    'rRNA', 'snoRNA', 'scaRNA', 'miRNA', 'snRNA', 'tRNA', 'Intron', 'sRNA',
    'CD-box', 'HACA-box',
    'CRISPR', 'antisense', 'leader', 'Cis-reg', 'Gene', 'splicing', ''
]
def process(bedfile, catalogue):
    familydescr = {}
    for line in open(catalogue):
        fields = line[:-1].split('\t')
        family = fields[0]
        descriptions = map(str.strip, fields[3].rstrip(';').split(';'))
        clsdescr = CLSPRIORITY[min(map(CLSPRIORITY.index, descriptions))]
        clsdescr = clsdescr if clsdescr != 'Gene' else 'ncRNA'
        familydescr[family] = clsdescr

    for line in bedfile:
        fields = line[:-1].split('\t')
        descrfields = fields[3].split(';')
        family = familydescr[fields[3].split(';')[0]]
        if descrfields[1].startswith('mir-'):
            continue
        fields[3] = '%s|%s|%s' % (family, descrfields[1], descrfields[1])
        print '\t'.join(fields[:6])

if __name__ == '__main__':
    import sys
    bedfile = sys.stdin
    catalogue = sys.argv[1]
#    bedfile = 'mm9/rnagene.bed'
#    catalogue = 'Rfam.catalogue'

    process(bedfile, catalogue)
