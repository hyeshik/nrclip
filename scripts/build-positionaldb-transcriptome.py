#!/usr/bin/env python
import bsddb
import numpy as np
import lzo as comp
import os

DEPTH_WINDOW = 10
THREEP_EXTENSION = 78
ROWSPACE = ['depth', '5', '3', 'A', 'C', 'G', 'T', 'I', 'D', 'N']
ROWSPACE2i = dict((n, i) for i, n in enumerate(ROWSPACE))
NROWS = len(ROWSPACE)
ARRSAVING = NROWS-1 # ignore N

def process(inpf, refflatdb, outdb):
    arrtracc = None
    arr = None

    def flush():
        if arrtracc is not None:
            outdb[arrtracc] = comp.compress(arr[:ARRSAVING].tostring())

    for line in inpf:
        fields = line[2:-1].split()
        nread = int(fields[0])
        tracc = fields[1]
        assert fields[2] == '+'
        pos = int(fields[3])
        mtype = line[0]

        if arrtracc != tracc:
            flush()
            arrtracc = tracc
            arrlength = refflatdb[tracc]['totalLength']
            arr = np.zeros([NROWS, arrlength + THREEP_EXTENSION], 'I')

        if mtype == 'P':
            tagspanleft = max(pos - DEPTH_WINDOW, 0)
            tagspanright = min(int(fields[4]) + DEPTH_WINDOW, arr.shape[1])
            arr[0, tagspanleft:tagspanright] += nread
        elif mtype == 'M':
            for relp, base in enumerate(fields[4]):
                p = pos + relp
                arr[ROWSPACE2i[base], p] += nread
        elif mtype in '53ID':
            arr[ROWSPACE2i[mtype], pos] += nread
        else:
            raise ValueError("Unknown message type %s" % repr(mtype))

    flush()

if __name__ == '__main__':
    import sys
    import shelve
    import os

    dboutpath = sys.argv[2]
    if os.path.exists(dboutpath):
        os.unlink(dboutpath)

    refflatdb = shelve.open(sys.argv[1], 'r')
    dbout = bsddb.hashopen(dboutpath, 'c')

    process(sys.stdin, refflatdb, dbout)

