#!/usr/bin/env python
TRANSCRIPTOME_MODE = False

def parse_chunk(chunk):
    chunkiter = iter(chunk)
    query = chunkiter.next()[1:-1].split('\t')

    targets = []
    for line in chunkiter:
        fields = line[1:-1].split('\t')
        strand = fields[2][0]
        if TRANSCRIPTOME_MODE and strand == '-':
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


def process(inf, outf):
    for orig, query, targets in bracket_chunk(inf):
        if len(targets) < 1:
            continue
        elif len(targets) == 1:
            best = [0]
        else:
            scores = [t['score'] for t in targets]
            minscore = min(scores)
            best = [i for i, s in enumerate(scores) if minscore == s]

        if len(best) == 1:
            print ''.join(orig)

if __name__ == '__main__':
    import sys

    if '-tr' in sys.argv:
        TRANSCRIPTOME_MODE = True
    process(sys.stdin, sys.stdout)
    #import gzip
    #process(gzip.open('../alignments/PolyA-Lin28a-clean.gmap.gz'), sys.stdout)
