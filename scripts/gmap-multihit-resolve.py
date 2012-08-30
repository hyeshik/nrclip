#!/usr/bin/env python
TRANSCRIPTOME_MODE = False
MAXIMUM_INDEL = 999999

def parse_aln_count(alns):
    rcnt = {}
    for chunk in [aln[3].split('..')[1] for aln in alns]:
        for token in chunk.split(','):
            key, cnt = token.split(':')
            if key not in rcnt:
                rcnt[key] = 0
            rcnt[key] += int(cnt)
    rcnt['indel'] = rcnt.get('ins', 0) + rcnt.get('del', 0)
    return rcnt

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

    filteredtargets = [
        tgt for tgt in targets
        if parse_aln_count(tgt['alns'])['indel'] <= MAXIMUM_INDEL
    ]
    #if len(targets) != len(filteredtargets):
    #    print '================='
    #    print targets
    #    print filteredtargets

    return query, filteredtargets

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
        sys.argv.remove('-tr')

    if len(sys.argv) >= 2:
        MAXIMUM_INDEL = int(sys.argv[1])

    process(sys.stdin, sys.stdout)
    #import gzip
    #process(gzip.open('../alignments/PolyA-Lin28a-clean.gmap.gz'), sys.stdout)
