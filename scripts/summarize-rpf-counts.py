#!/usr/bin/env python
from __future__ import division
import sys
import csv
from itertools import chain

TOTALREADS = ''
NORMBASE = 1000000.

def load_counts(filename):
    r = {}
    for line in open(filename):
        name, cnt = line.split()
        r[name] = int(cnt)

    r[TOTALREADS] = sum(r.itervalues())
    return r

def process(readcounts, samplepairs):
    alltranscripts = reduce(lambda x, y: x | set(y),
                            readcounts.itervalues(), set())
    w = csv.writer(sys.stdout, dialect='excel-tab')
    w.writerow(['accession'] +
               list(chain(*[['%s PolyA' % sample, '%s RPF' % sample,
                             '%s density' % sample]
                            for sample, _, _ in samplepairs])))

    def normalized_count(label, accession):
        return ((readcounts[label].get(accession, 0) + 1) * NORMBASE /
                readcounts[label][TOTALREADS])
        #return readcounts[label].get(accession, 0)

    for acc in sorted(alltranscripts - set([TOTALREADS])):
        fields = [acc]
        for sample, polya_label, rpf_label in samplepairs:
            polya_value = normalized_count(polya_label, acc)
            rpf_value = normalized_count(rpf_label, acc)
            fields.extend([polya_value, rpf_value, rpf_value / polya_value])
            #fields.extend([polya_value, rpf_value, (rpf_value+1) / (polya_value+1)])
        w.writerow(fields)


if __name__ == '__main__':
#    cntfilemapping = """
#        1R:../stats/cdsreads.RPF-Luc-rep1.txt
#        1P:../stats/cdsreads.PolyA-Luc-rep1.txt
#        2R:../stats/cdsreads.RPF-Luc-rep2.txt
#        2P:../stats/cdsreads.PolyA-Luc-rep2.txt
#    """
#    samplepairs = "Luc-r1:1P:1R Luc-r2:2P:2R"
    cntfilemapping = sys.argv[1]
    samplepairs = sys.argv[2]

    cntfilemapping = dict(token.split(':') for token in cntfilemapping.split())
    samplepairs = [
        tuple(token.split(':')) for token in samplepairs.split()]

    readcounts = dict((label, load_counts(filename))
                      for label, filename in cntfilemapping.iteritems())

    process(readcounts, samplepairs)

