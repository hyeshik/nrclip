#!/usr/bin/env python
from __future__ import division
import numpy as np
import pickle

def calc_statistics(table, width):
    epartedr = table.transpose().reshape([width, 5, 5])

    positional_total_reads = []
    for i in range(width):
        positional_total_reads.append(epartedr[i].sum(axis=0))

    positional_total_reads = np.array(positional_total_reads)

    proportions = np.zeros(epartedr.shape, 'f8')
    for i in range(width):
        #print '-----', i
        proportions[i] = epartedr[i] / positional_total_reads[i].clip(1)
        #print np.round(proportions[i] * 100, 1)

    subst_rate = np.zeros([width, 4], 'f8')
    del_rate = np.zeros([width, 4], 'f8')
    ins_rate = np.zeros([width, 4], 'f8')
    for i in range(width):
        for j, nt in enumerate('ACGT'):
            subst_rate[i, j] = proportions[i, :4, j].sum(0) - proportions[i, j, j]
            del_rate[i, j] = proportions[i, 4, j]

        ins_rate[i] = epartedr[i, :4, 4] / epartedr[i].sum(axis=1)[:4].clip(1)

    matchmask = np.identity(5, dtype='i8')
    substmask = np.ones(matchmask.shape, dtype='i8') - matchmask
    delmask = substmask.copy()
    insmask = substmask.copy()

    matchmask[:, 4] = substmask[:, 4] = substmask[4, :] = \
        delmask[:4, :] = insmask[:, :4] = 0

    any_subst_rate = []
    any_del_rate = []
    any_ins_rate = []

    for i in range(width):
        nmatch = (epartedr[i, :] * matchmask).sum()
        nsubst = (epartedr[i, :] * substmask).sum()
        ndel = (epartedr[i, :] * delmask).sum()
        nins = (epartedr[i, :] * insmask).sum()

        total = nmatch + nsubst + ndel + nins
        any_subst_rate.append(nsubst / total)
        any_del_rate.append(ndel / total)
        any_ins_rate.append(nins / total)

    return {
        'subst': subst_rate,
        'del': del_rate,
        'ins': ins_rate,
        'anysubst': any_subst_rate,
        'anydel': any_del_rate,
        'anyins': any_ins_rate,
    }

def process(inpickle):
    eleft, eright, eparted = pickle.load(inpickle)

    left = calc_statistics(eleft[:, 0:78], 78)
    right = calc_statistics(eright[:, 0:78], 78)
    parted = calc_statistics(eparted, 20)
    return left, right, parted

if __name__ == '__main__':
    import sys
    import pickle

    readerrorfile = sys.argv[1]
    outfile = sys.argv[2]

    r = process(open(readerrorfile))
    pickle.dump(r, open(outfile, 'w'))

