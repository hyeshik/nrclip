#!/usr/bin/env python
from __future__ import division
import gzip
import os
import numpy as np
import re
from itertools import groupby
from nrclipext import shannon_entropy
from rnarry.sequtils import GiantFASTAFile

BASE2i = dict((b, i) for i, b in enumerate('ACGT-'))
LONGDELETION = re.compile('-{3,}')

def parse_chunk(chunk):
    chunkiter = iter(chunk)
    query = chunkiter.next()[1:-1].split('\t')

    targets = []
    for line in chunkiter:
        fields = line[1:-1].split('\t')
        strand = fields[2][0]
        if strand == '-':
            pass
        elif line[0] == ' ':
            score = [int(token[12:]) for token in fields[4].split(',')
                     if token.startswith('align_score:')][0]
            targets.append({'score': score, 'alns': [fields]})
        elif line[0] == ',':
            targets[-1]['alns'].append(fields)

    return query, targets

def bracket_chunk(inf):
    r = []
    for line in inf:
        if line.startswith('>'):
            if r:
                yield (r,) + parse_chunk(r)
            r = [line]
        elif line != '\n':
            r.append(line)
    if r:
        yield (r,) + parse_chunk(r)

def scan_gmap(mapinp):
    spans = []

    for orig, query, targets in bracket_chunk(mapinp):
        nreads = int(query[2].split('-')[1])
        qflen = len(query[0])
        fwd_targets = [tgt for tgt in targets if tgt['alns'][0][2][0] == '+']
        if not fwd_targets:
            continue

        alns = fwd_targets[0]['alns']
        tname = alns[0][2][1:].split(':')[0]
        tstarts = [int(aln[2].split(':')[1].split('..')[0]) - 1
                   for aln in alns]

        if len(alns) == 1 and alns[0][0] == query[0]:
            # fast path for perfect matches
            spans.append((tname, tstarts[0], nreads, query[0]))
            continue

        tends = [int(aln[2].split(':')[1].split('..')[1]) for aln in alns]
        qstarts = [int(aln[1].split('..')[0]) - 1 for aln in alns]
        qends = [int(aln[1].split('..')[1]) for aln in alns]

        seqassigned = {}

        tpositions = tstarts[:]
        segseqs = [a[0] for a in alns]
        targetalnrange = range(1, len(alns)+1)

        for i, cluster in enumerate(zip(*[query[0]] + segseqs)):
            refread = cluster[0]
            for j in range(1, len(alns)+1):
                r = cluster[j]
                if r == '-':
                    continue

                talni = j - 1
                if qstarts[talni] <= i < qends[talni]:
                    tpos = tpositions[talni]
                    if tstarts[talni] <= tpos < tends[talni]:
                        seqassigned[tpos] = refread.upper()
                    tpositions[talni] += 1

        reconstitutedseq = ''.join(seqassigned.get(i, '-')
                                    for i in range(tstarts[0], tends[-1]))
        if not LONGDELETION.findall(reconstitutedseq):
            spans.append((tname, tstarts[0], nreads, reconstitutedseq))

    spans.sort()
    return spans

def probe_reads(spans, greffile, nonzeroout, distout):
    outd = distout['D']
    outm = distout['M']
    outmd = distout['MD']
    oute = distout['E']

    for spname, grp in groupby(spans, key=lambda sp: sp[0]):
        seq = greffile.get(spname).upper()
        posreads = np.zeros([len(seq), 5], 'I')

        for _, start, nreads, readseq in grp:
            for offset, r in enumerate(readseq):
                if r in BASE2i:
                    posreads[start+offset, BASE2i[r]] += nreads

        nonzeropositions = np.where(posreads.sum(1) > 0)[0]
        for pos in nonzeropositions:
            refbase = seq[pos]
            simcnt = posreads[pos]
            kcnt = simcnt.sum()
            del_ratio = simcnt[4] / kcnt
            if refbase in BASE2i:
                refcnt = simcnt[BASE2i[refbase]]
                moddel_ratio = (kcnt - refcnt) / kcnt
                mod_ratio = (kcnt - refcnt - simcnt[4]) / kcnt
            else:
                mod_ratio = moddel_ratio = 0.

            sentropy = shannon_entropy(simcnt)

            if 0 < del_ratio < 0.9 or 0 < mod_ratio < 0.9:
                print >> nonzeroout, spname, pos, refbase, kcnt, ' '.join(
                 map(str, simcnt)), del_ratio, mod_ratio, moddel_ratio, sentropy

            print >> outd, kcnt, refbase, '%.6f' % del_ratio
            print >> outm, kcnt, refbase, '%.6f' % mod_ratio
            print >> outmd, kcnt, refbase, '%.6f' % moddel_ratio
            print >> oute, kcnt, refbase, '%.6f' % sentropy


if __name__ == '__main__':
    import sys

    mappingfile = sys.argv[1] # gmap native format
    greffile = sys.argv[2]
    outputprefix = sys.argv[3]

    mappingf = gzip.open(mappingfile)
    greffile = GiantFASTAFile(greffile)

    def opensortedout(fname):
        return os.popen("sort -k1,1n -k2,2 -k3,3nr | uniq -c | "
                "awk '{ print $2, $3, $4, $1; }' | gzip - >%s" % fname, 'w')

    nonzeroout = gzip.open(outputprefix + 'nonzero.posrcnt.gz', 'w')
    distout = {
        'D': opensortedout(outputprefix + 'del.real.gz'),
        'M': opensortedout(outputprefix + 'mod.real.gz'),
        'MD': opensortedout(outputprefix + 'moddel.real.gz'),
        'E': opensortedout(outputprefix + 'entropy.real.gz'),
    }

    spans = scan_gmap(mappingf)
    probe_reads(spans, greffile, nonzeroout, distout)

