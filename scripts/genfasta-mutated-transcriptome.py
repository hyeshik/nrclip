#!/usr/bin/env python
import shelve
import numpy as np
import sys
from rnarry.sequtils import GiantFASTAFile
from rnarry.sequtils import reverse_complement
from rnarry.utils import textwrap
from rnarry.bxwrap import MultiTrackSplitBinnedArray

MINREADSTOCALL = 10
MINPERCENTTOCALL = 0.9

blklength = lambda blks: sum(end-start for start, end in blks)

class ContinuousCache(object):

    def __init__(self, opener):
        self.curkey = None
        self.current = None
        self.opener = opener

    def get(self, key):
        if self.curkey != key:
            self.curkey = key
            self.current = self.opener(key)

        return self.current

def process(nrlist):
    nucleotidetracks = [rseqarr.TRACKS.index(i) for i in 'ACGT']

    for nracc in nrlist:
        dbinfo = refFlat[nracc]
        chrom = dbinfo['chrom']
        genomeseq = ''.join(mm9.get(chrom, blkstart, blkend)
                            for blkstart, blkend in dbinfo['exonBlocks']).upper()
        if dbinfo['strand'] == '-':
            utr3, utr5 = 'leftUtrBlocks', 'rightUtrBlocks'
        else:
            utr5, utr3 = 'leftUtrBlocks', 'rightUtrBlocks'

        if nracc.startswith('NM_'):
            utr5length = blklength(dbinfo[utr5])
            cdslength = blklength(dbinfo['cdsBlocks'])
            utr3length = blklength(dbinfo[utr3])
        else:
            exonlength = blklength(dbinfo['exonBlocks'])

        cntarray = rseqarr.get_blocks(dbinfo['chrom'], dbinfo['exonBlocks'],
                                      dbinfo['strand'])[nucleotidetracks]
        depthcnt = np.array(cntarray.sum(0).clip(1), 'd')
        confidentcalls = ((cntarray/depthcnt >= MINPERCENTTOCALL) *
                          (depthcnt >= MINREADSTOCALL))

        mutatedseq = list(genomeseq)
        for base, calls in zip('ACGT', confidentcalls):
            for pos in np.where(calls)[0]:
                mutatedseq[pos] = base

        mutatedseq = ''.join(mutatedseq)

        if dbinfo['strand'] == '-':
            mutatedseq = reverse_complement(mutatedseq)

        if nracc.startswith('NM_'):
            print >> bedout, '\t'.join([nracc, str(utr5length),
                                        str(utr5length + cdslength),
                                        '%s' % dbinfo['geneName'], '.', '+'])
        else:
            print >> bedout, '\t'.join([nracc, '0', str(exonlength),
                                        '%s' % dbinfo['geneName'], '.', '+'])

        print >> faout, '>%s %s' % (nracc, dbinfo['geneName'])
        faout.write(textwrap(mutatedseq))


if __name__ == '__main__':
    refflatdbpath = sys.argv[1]
    nrlistpath = sys.argv[2]
    genomefastapath = sys.argv[3]
    rnaseqgspace = sys.argv[4]
    fastaoutpath = sys.argv[5]
    bedoutpath = sys.argv[6]

    refFlat = shelve.open(refflatdbpath, 'r')
    mm9 = GiantFASTAFile(genomefastapath)
    rseqarr = MultiTrackSplitBinnedArray(rnaseqgspace)
    nrlist = open(nrlistpath).read().split()

    bedout = open(bedoutpath, 'w')
    faout = open(fastaoutpath, 'w')

    process(nrlist)

