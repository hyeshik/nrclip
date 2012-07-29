#!/usr/bin/env python
from __future__ import division
import cPickle as pickle
import numpy as np
import math
import shelve
from itertools import chain
from rpy import r
import csv

pchisq = r.pchisq
padjust = r.p_adjust
RATIO_MINREADS = 50
SPLITCNTLABELS = ('5UTR', 'CDS', '3UTR', 'Exon')
NULLREADS = (0, 0, 0)

def load(d, name, infile):
    splitcnts = d
    splitcntidx = dict((label, i) for i, label in enumerate(SPLITCNTLABELS))

    p = pickle.load(open(infile))
    for (acc, label), cnts in p.iteritems():
        cntplace = splitcnts[splitcntidx[label]]
        if acc not in cntplace:
            cntplace[acc] = {}
        cntplace[acc][name] = cnts

    return splitcnts

def calclogratio(d, samplea, sampleb):
    r = {}
    for acc, cnts in d.iteritems():
        cnta = cnts.get(samplea, NULLREADS)
        cntb = cnts.get(sampleb, NULLREADS)
        if cnta[0] < RATIO_MINREADS and cntb[0] < RATIO_MINREADS:
            continue

        l2rat = math.log((cntb[0] + 1) / (cnta[0] + 1), 2)
        r[acc] = l2rat

    return r

def writeresult(refsample, samples, origcnt, l2ratio, significance):
    w = csv.writer(sys.stdout)
    w.writerow(['accession', 'gene', 'p-value', 'FDR',
                'exonCount', 'totalLength', 'cdsLength', 'noncodingLength'] +
                list(chain(*[['%s 5UTR cnt' % sample, '%s CDS cnt' % sample,
                              '%s 3UTR cnt' % sample, '%s Exon cnt' % sample,
                              '%s 5UTR l2r' % sample, '%s CDS l2r' % sample,
                              '%s 3UTR l2r' % sample, '%s Exon l2r' % sample]
                              for sample in samples])))

    l2idx = range(len(SPLITCNTLABELS))

    allacc = reduce(lambda x, y: x|y, (set(cnt) for cnt in origcnt))
    for acc in allacc:
        if all((acc not in inst)
                for sample in samples for inst in l2ratio[sample]):
            continue

        entry = REFFLAT[acc]
        if 'cdsBlocks' in entry:
            cdslen = sum(b-a for a, b in entry['cdsBlocks'])
        else:
            cdslen = 0

        gsym = entry['geneName']
        out = [acc, gsym, significance.get(acc, (1., 1.))[0],
               significance.get(acc, (1., 1.))[1], entry['exonCount'],
               entry['totalLength'], cdslen, entry['totalLength'] - cdslen]
        for sample in samples:
            out.extend([cnt.get(acc, {}).get(sample, NULLREADS)[0]
                        for cnt in origcnt])
            out.extend([rat.get(acc, 0.) for rat in l2ratio[sample]])

        w.writerow(out)

def calc_significance(ratios):
    # ranking for nonparametric tests
    commonaccs = reduce(lambda x, y: x & y,
                        (set(ratios[smp][3]) for smp in SAMPLES_FOR_FDR))
    accranks = None
    for sample in SAMPLES_FOR_FDR:
        sratios = ratios[sample]
        ordered = sorted(commonaccs,
                         key=lambda acc: sratios[3][acc], reverse=True)
        if accranks is None:
            accranks = dict((acc, [rank]) for rank, acc in enumerate(ordered))
        else:
            for rank, acc in enumerate(ordered):
                accranks[acc].append(rank)

    # chi-squared test
    totaltranscripts = len(accranks)

    accs, pvalues = [], []
    for acc, ranks in accranks.iteritems():
        chisq = -2 * np.log(np.add(ranks, 1) / totaltranscripts).sum()
        pvalue = 1 - pchisq(chisq, len(SAMPLES_FOR_FDR)*2)
        accs.append(acc)
        pvalues.append(pvalue)

    fdrs = padjust(pvalues, 'BH')

    return dict(zip(accs, zip(pvalues, fdrs)))

if __name__ == '__main__':
    import sys

    REFFLAT = shelve.open(sys.argv[1], 'r')
    refsample = sys.argv[2]
    clipsamples = sys.argv[3].split(',')
    datafiles = [token.split(':') for token in sys.argv[4:]]
    samples = [label for label, _ in datafiles]

    SAMPLES_FOR_FDR = clipsamples

    d = {}, {}, {}, {}

    for label, fname in datafiles:
        load(d, label, fname)

    # add 0, 1, 2 (UTRs and CDS) to 3 (exon)
    exoncnt = d[3]
    for partcnt in d[0:3]:
        for acc, smpcnt in partcnt.iteritems():
            if acc not in exoncnt:
                exoncnt[acc] = smpcnt
                continue

            currcnt = exoncnt[acc]
            for smp, cnt in smpcnt.iteritems():
                currcnt[smp] = tuple(np.add(cnt, currcnt[smp]))

    ratios = {}
    for sample in samples:
        l2r = []
        for iidx in range(len(d)):
            l2r.append(calclogratio(d[iidx], refsample, sample))
        ratios[sample] = tuple(l2r)

    significance = calc_significance(ratios)

    writeresult(refsample, samples, d, ratios, significance)

