#!/usr/bin/env python
import numpy as np

MAXREADLENGTH = 100
LABELS = """
AA CA GA TA IA
AC CC GC TC IC
AG CG GG TG IG
AT CT GT TT IT
AD CD GD TD XX""".split()
LABEL2i = dict((l, i) for i, l in enumerate(LABELS))

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


def process(inf, readlengthlim):
    matchcount = np.zeros([len(LABELS), MAXREADLENGTH], 'i8')
    matchcount_r = np.zeros([len(LABELS), MAXREADLENGTH], 'i8')
    matchcount_20 = np.zeros([len(LABELS), 20], 'i8')

    def addmatch(label, i, l, nr, ign=set(['NI', 'IN', 'AN', 'CN', 'TN', 'GN',
                                        'NG', 'NC', 'NT', 'NA'])):
        #print ' -- <%s:%d> += %d' % (label, i, nr)
        if label not in ign:
            li = LABEL2i[label]
            matchcount[li, i] += nr
            if l <= readlengthlim:
                matchcount_r[li, max(0, l-i-1)] += nr
                matchcount_20[li, i * 20 // l] += nr

    for orig, query, targets in bracket_chunk(inf):
        nreads = int(query[2].split('-')[1])
        qflen = len(query[0])
        #print query, nreads, targets[0]['alns'][0][-2]

        alns = targets[0]['alns']
        if len(alns) == 1 and alns[0][0] == query[0]:
            # fast path for perfect matches
            for i, nt in enumerate(query[0]):
                addmatch(nt*2, i, qflen, nreads)
            continue

        segseqs = [a[0] for a in alns]
        lastmatches = [999] * len(segseqs)
        for i, cluster in enumerate(zip(*[query[0]] + segseqs)):
            nonempty = filter(lambda x: x!= '-', cluster)
            for j, tg in enumerate(cluster[1:]):
                if tg.isupper(): # match
                    addmatch(tg*2, i, qflen, nreads)
                    lastmatches[j] = i
                elif tg.islower():
                    if len(nonempty) <= 2: # substitution
                        addmatch(tg.upper() + cluster[0], i, qflen, nreads)
                        lastmatches[j] = i
                    else: # deletion
                        addmatch(tg.upper() + 'D',
                                 min(i, lastmatches[j]+1), qflen, nreads)
                elif tg == '-' and len(nonempty) == 1: # insertion
                    addmatch('I' + cluster[0], i, qflen, nreads)
                    break

    return matchcount, matchcount_r, matchcount_20

if __name__ == '__main__':
    import sys
    import gzip
    import pickle

    gmapfile = sys.argv[1]
    outputfile = sys.argv[2]
    readlengthlim = int(sys.argv[3])

    r = process(gzip.open(gmapfile), readlengthlim)
    pickle.dump(r, open(outputfile, 'w'))

