#!/usr/bin/env python
#
# Copyright (c) 2012 Hyeshik Chang.
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# - Hyeshik Chang <hyeshik@snu.ac.kr>
#

from Bio import SeqIO
import sys
from collections import defaultdict
import csv

def collect_statistics(inputfile):
    totalcount = defaultdict(int)
    ntcount = dict((nt, defaultdict(int)) for nt in 'ACGTN')

    for seq in SeqIO.parse(inputfile, format='fasta'):
        nreads = int(seq.description.split('-')[1])
        totalcount[len(seq)] += nreads
        if len(seq) > 0:
            ntcount[str(seq.seq)[-1]][len(seq)] += nreads

    w = csv.writer(sys.stdout)
    w.writerow(['length', 'reads', '-A', '-C', '-G', '-T', '-N'])

    for i in sorted(totalcount):
        w.writerow([i, totalcount[i]] +
                   [ntcount[nt][i] for nt in 'ACGTN'])

if __name__ == '__main__':
    filename = sys.argv[1]
    collect_statistics(open(filename))

