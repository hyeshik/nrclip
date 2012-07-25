#!/usr/bin/env python
from bx.binned_array import BinnedArray
import os

COMPTYPE = 'zlib'
TYPECODE = 'I'
newarray = lambda: BinnedArray(default=0, typecode=TYPECODE)

def open_with_newdir(filename):
    dname = os.path.dirname(filename)
    if not os.path.isdir(dname):
        os.makedirs(dname)
    return open(filename, 'w')

BASESPACE = ['depth', 'A', 'C', 'G', 'T', 'N']
BASESPACE2i = dict((n, i) for i, n in enumerate(BASESPACE))

def process(inpf, prefix):
    arrchrom = None
    arrtype = None
    arr = []

    def flush():
        if arrchrom is None:
            return

        filenamefmt = '%s/%%s/%s.%s' % (prefix, arrchrom[0], arrchrom[1])
        if arrtype != 'M':
            arr[0].to_file(open_with_newdir(filenamefmt % arrtype), COMPTYPE)
        else:
            for subtype, a in zip(BASESPACE, arr):
                a.to_file(open_with_newdir(filenamefmt % subtype), COMPTYPE)
        del arr[:]

    for line in inpf:
        fields = line[2:-1].split()
        nread = int(fields[0])
        chrom = fields[1], fields[2]
        pos = int(fields[3])
        mtype = line[0]

        if arrchrom != chrom or arrtype != mtype:
            flush()
            arrchrom = chrom
            arrtype = mtype
            if mtype != 'M':
                arr = [newarray()]
            else:
                arr = [newarray() for i in BASESPACE]

        if mtype == 'M':
            for relp, base in enumerate(fields[4]):
                p = pos + relp
                barr = arr[BASESPACE2i[base]]
                arr[0].set(p, arr[0].get(p) + nread)
                barr.set(p, barr.get(p) + nread)
        elif mtype in '53ID':
            arr[0].set(pos, arr[0].get(pos) + nread)
        else:
            raise ValueError("Unknown message type %s" % repr(mtype))

    flush()

    open('%s/done' % prefix, 'w')

if __name__ == '__main__':
    import sys
    process(sys.stdin, sys.argv[1])

