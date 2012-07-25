#
# rnarry.nrclip.ContaminantFilter
#  - Removes contaminants like rRNA and tRNA
#
#
# Copyright (C) 2012 Hyeshik Chang
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
#

from ruffus import *
from rnarry.nrclip import Paths, Options, SequenceProcessing, DataPreparation
from rnarry.nrclip.PipelineControl import *


@files(for_each(Paths.fulltag_collapsed_reads,
                Paths.fulltag_prealn_sam,
                Paths.ALL_SAMPLES,
                [Options.FULLTAG_PREALN_MISMATCHES]) +
       for_each(Paths.shorttag_tags,
                Paths.shorttag_prealn_sam,
                Paths.SHORTTAG_SAMPLES,
                [Options.SHORTTAG_PREALN_MISMATCHES]))
@follows(SequenceProcessing.fulltag_collapse)
@follows(SequenceProcessing.shorttag_trim_collapse)
@follows(DataPreparation.generate_gsnap_early_filter_index)
@jobs_limit(1, 'exclusive')
def align_to_known_contaminants(inputfile, outputfile, sample):
    runproc("""
        $GSNAP -D $EXTERNAL_DIR -d $early_filter_prefix -B 4 -A sam \
            -m $FULLTAG_PREALN_MISMATCHES -t $NUM_THREADS $inputfile | \
        $GZIP -c - > $outputfile""", outputfile)


@files(for_each((Paths.fulltag_collapsed_reads, Paths.fulltag_prealn_sam),
                Paths.fulltag_filtered_reads,
                Paths.ALL_SAMPLES) +
       for_each((Paths.shorttag_tags, Paths.shorttag_prealn_sam),
                Paths.shorttag_filtered_reads,
                Paths.SHORTTAG_SAMPLES))
@follows(align_to_known_contaminants)
def filter_contaminants(inputfiles, outputfile, sample):
    seq_fasta, prealn_sam = inputfile

    with TemporaryFile() as idlist:
        # prepare list of sequence IDs
        runproc("""
            $ZCAT $prealn_sam | \
            $GREP -v '^@' | \
            $AWK -F' ' '{ if ($$3 != "*") { print $$1; } }' | \
            $SORT -t- -k1,1n | $UNIQ > $idlist""")

        # filter sequences
        runproc('$FASOMERECORDS -exclude $seq_fasta $idlist $outputfile',
                outputfile)


def tasks():
    return [
        align_to_known_contaminants,
        filter_contaminants,
    ]

