#!/usr/bin/env python
import gzip
import re
from rnarry.interfaces import sam
from rnarry.sequtils import reverse_complement

ENDTRIM = 2

cigar_pattern = re.compile('(\d+)([MIDNSHP])')

def genrefblks(readseq, chrom, start, stop, strand, cigar, nreads):
    refpos = start
    readpos = 0
    if strand == '-':
        readseq = reverse_complement(readseq)

    tleftlim, trightlim = start + ENDTRIM, stop - ENDTRIM
    qleftlim, qrightlim = ENDTRIM, len(readseq) - ENDTRIM
    #print '==============='
    #print readseq, chrom, start, stop, strand, cigar, nreads

    if transcriptome:
        print 'P', nreads, chrom, strand, tleftlim, trightlim

    cigarcommands = cigar_pattern.findall(cigar)
    if cigarcommands[0][1] == 'S': # shift start site for the first soft clipping
        start -= int(cigarcommands[0][0])
    if len(cigarcommands) > 1 and cigarcommands[-1][1] == 'S': # last soft clipping
        stop += int(cigarcommands[-1][0])

    for num, cmd in cigarcommands:
        num = int(num)
        if cmd == 'M': # match
            mleft = max(qleftlim, readpos)
            mright = min(qrightlim, readpos + num)
            if mleft < mright:
                seq = readseq[mleft:mright]
                print 'M', nreads, chrom, strand, max(refpos, tleftlim), seq
            refpos += num
            readpos += num
        elif cmd == 'S': # soft clip
            readpos += num
        elif cmd == 'N': # skip
            refpos += num
        elif cmd == 'D': # deletion
            if tleftlim <= refpos < trightlim:
                print 'D', nreads, chrom, strand, refpos, num
            refpos += num
        elif cmd == 'I': # insertion
            ppos = (refpos if strand == '+' else (refpos-1))
            if tleftlim <= ppos < trightlim:
                print 'I', nreads, chrom, strand, ppos, num
            readpos += num
        elif cmd == 'H': # hard clipping
            pass
        else:
            print 'E', nreads, num, cmd, readseq
            raise ValueError

    if strand == '+':
        fivep, threep = start, stop-1
    else:
        fivep, threep = stop-1, start

    print '5', nreads, chrom, strand, fivep
    print '3', nreads, chrom, strand, threep

def process(insam):
    for aln in insam:
        nreads = int(aln['qname'].split('-')[1])
        for chrom, start, stop, strand, mm, cigar in aln['mapped']:
            genrefblks(aln['seq'], chrom, start, stop, strand, cigar, nreads)

if __name__ == '__main__':
    import sys
    inpf = gzip.open(sys.argv[1])
    if len(sys.argv) > 2 and sys.argv[2] == '-tr':
        transcriptome = True
    else:
        transcriptome = False

    insam = sam.SAMParser(inpf, zerobase=True)
    process(insam)

