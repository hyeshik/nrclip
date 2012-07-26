#!/usr/bin/env python
import cPickle as pickle
from itertools import groupby
import sys
import os

MINDEPTH = int(os.environ['MINDEPTH'])

def load(expgenes, picklef):
    cntdata = pickle.load(picklef)
    for acc, grp in groupby(sorted(cntdata.items()), key=lambda x: x[0][0]):
        totalread = sum(int(cnt[0]) for k, cnt in grp)
        if totalread < MINDEPTH:
            continue

        expgenes.add(acc)

if __name__ == '__main__':
    nrlist = set(open(sys.argv[1]).read().split())
    expgenes = set()
    for f in sys.argv[2:]:
        load(expgenes, open(f))
    print '\n'.join(sorted(expgenes & nrlist))

