#!/usr/bin/env python
from __future__ import division
import pickle
import gzip
import os
import struct
import glob
import numpy as np
import re
from itertools import groupby
from rnarry.sequtils import GiantFASTAFile
from rnarryext import RandomReadCounter, shannon_entropy

def detect_readcount_type():
    if struct.calcsize('I') != 4:
        raise SystemError("`int' type needs to be a 32-bit primitive type.")

    for tp in 'LQ':
        if struct.calcsize(tp) == 8:
            return tp
    raise SystemError('No 64-bit types available in struct module')

READCOUNT_TYPE = detect_readcount_type()
i2BASE = 'ACGT-'
BASE2i = dict((b, i) for i, b in enumerate(i2BASE))

def load_profile(inpf):
    prof, _, _ = pickle.load(open(inpf))

    maxbasepos = prof.shape[1]
    positional_profile = prof.transpose().reshape([maxbasepos, 5, 5])
    maxbasepos = max(i for i in range(maxbasepos)
                     if positional_profile[i].sum() > 0)

    return positional_profile, maxbasepos


def dump_profile(inpprofile, maxpos, output):
    output.write(struct.pack('II', maxpos+1, len(i2BASE)))

    profile_struct = READCOUNT_TYPE * 5
    for pos in range(maxpos+1):
        for ni in range(len(i2BASE)):
            readdist = inpprofile[pos, :, ni]
            output.write(struct.pack(profile_struct, *list(readdist)))
            #if pos == 0:
            #    print readdist


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
        spanname = spanstart = spanend = None
        nreads = int(query[2].split('-')[1])
        for aln in targets[0]['alns']:
            if aln[2][0] == '-':
                continue
            tname, trange = aln[2][1:].split(':')
            tstart, tend = map(int, trange.split('..'))
            tstart -= 1

            if spanstart is None:
                spanstart = tstart
            spanend = tend

            assert spanname is None or spanname == tname
            spanname = tname

        spans.append((spanname, spanstart, spanend, nreads))

    spans.sort()
    return spans


def dump_gmap(spans, output, greffile):
    for name, grp in groupby(spans, lambda x: x[0]):
        grp = list(grp)
        if len(name) > 31:
            raise ValueError("transcript name %s too long." % repr(name))

        refseq = greffile.get(name)
        output.write(name.ljust(32, '\x00'))
        output.write(struct.pack('II', len(refseq), len(grp)))
        output.write(refseq)
        for _, start, end, nreads in grp:
            output.write(struct.pack('III', start, end, nreads))


if __name__ == '__main__':
    import sys

    refseqfile = sys.argv[1] # fasta format
    mappingfile = sys.argv[2] # gzipped gmap native format
    readerrorfile = sys.argv[3] # python pickle format
    packfile = sys.argv[4]

    greffile = GiantFASTAFile(refseqfile)
    mappingf = gzip.open(mappingfile)

    packf = open(packfile, 'w')

    print 'Loading read error profile'
    profile, maxpos = load_profile(readerrorfile)

    print 'Dumping read error profile'
    dump_profile(profile, maxpos, packf)

    print 'Loading gmap mapping file'
    spans = scan_gmap(mappingf)

    print 'Dumping mapped reads'
    dump_gmap(spans, packf, greffile)

