#
# rnarry.nrclip.SequenceProcessing
#  - Basic sequence processing routines
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
from rnarry.nrclip import Paths, Options
from rnarry.nrclip.PipelineControl import *

FASTX_TOOLKIT_QUALITY = {
    'illumina1.5': '',
    'illumina1.8': '-Q33',
}


@files(for_each(Paths.original_sequence_reads,
                Paths.fulltag_quality_filtered_reads,
                Paths.ALL_SAMPLES))
def fulltag_filter_clip_trim(inputfile, outputfile, sample):
    adapter = Options.ADAPTER_SEQ[sample]
    qoption = FASTX_TOOLKIT_QUALITY[Options.QUALITY_SCALE[sample]]
    runproc("""
        $ZCAT $inputfile | \
        $FASTX_CLIPPER %s -n -a $adapter -l $FULLTAG_MIN_LENGTH | \
        $FASTQ_QUALITY_TRIMMER %s -t $FULLTAG_MIN_QUALITY \
            -l $FULLTAG_MIN_LENGTH | \
        $FASTQ_QUALITY_FILTER %s -q $FULLTAG_MIN_QUALITY \
            -p $FULLTAG_MIN_QUALITY_PERCENT -z -o $outputfile"""
            % (qoption, qoption, qoption),
        outputfile)


@files(for_each(Paths.fulltag_quality_filtered_reads,
                Paths.fulltag_collapsed_reads,
                Paths.ALL_SAMPLES))
@follows(fulltag_filter_clip_trim)
def fulltag_collapse(inputfile, outputfile, sample):
    qoption = FASTX_TOOLKIT_QUALITY[Options.QUALITY_SCALE[sample]]
    runproc("""
        $ZCAT $inputfile | \
        $FASTX_COLLAPSER %s | \
        $FASTX_ARTIFACTS_FILTER -o $outputfile""" % qoption, outputfile)


@files(for_each(Paths.original_sequence_reads,
                Paths.shorttag_tags,
                Paths.SHORTTAG_SAMPLES))
def shorttag_trim_collapse(inputfile, outputfile, sample):
    qoption = FASTX_TOOLKIT_QUALITY[Options.QUALITY_SCALE[sample]]
    runproc("""
        $ZCAT $inputfile | \
        $FASTX_TRIMMER %s -f 1 -l $SHORTTAG_LENGTH | \
        $FASTQ_QUALITY_FILTER %s -q $SHORTTAG_MIN_QUALITY \
            -p $SHORTTAG_MIN_QUALITY_PERCENT | \
        $FASTX_COLLAPSER %s | \
        $FASTX_ARTIFACTS_FILTER -o $outputfile"""
            % (qoption, qoption, qoption), outputfile)


def tasks():
    return [
        fulltag_filter_clip_trim,
        fulltag_collapse,
        shorttag_trim_collapse,
    ]

