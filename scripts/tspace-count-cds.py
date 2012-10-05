#!/usr/bin/env python
from __future__ import division
import sys
import csv
from rnarry.intervalarray import TranscriptLevelReadsStats
import shelve

FOOTPRINTSHIFT = -15

NRREFSEQ_DB = sys.argv[1]
TSPACE_INPUT = sys.argv[2]

nrRefSeq = shelve.open(NRREFSEQ_DB, 'r')
tspace = TranscriptLevelReadsStats(TSPACE_INPUT)

cdsreads = {}

for acc in tspace.keys():
    if acc.startswith('NR_') or acc not in nrRefSeq:
        continue

    # get 5'-end counts array
    fpendscnt = tspace.get(acc, '5')

    # calculate cds interval
    refinfo = nrRefSeq[acc]
    utr5, cds, utr3 = refinfo['partLengths']
    orfstart = max(0, utr5 + FOOTPRINTSHIFT)
    orfend = min(utr5 + cds + FOOTPRINTSHIFT, utr5 + cds + utr3)

    cdsreads[acc] = int(fpendscnt.sum())

cdsreads = cdsreads.items()
cdsreads.sort(key=lambda (k, v): (-v, k))

w = csv.writer(sys.stdout, dialect='excel-tab')
w.writerows(cdsreads)
