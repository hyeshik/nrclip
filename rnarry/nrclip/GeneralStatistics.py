#
# rnarry.nrclip.GeneralStatistics
#  - Produces general statistics for sequencing quality control and overview
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
import os
from rnarry.nrclip import (
    Paths, Options, DerivedDatabaseBuilding, SequenceAnnotation,
    ContaminantFilter, SequenceProcessing)
from rnarry.nrclip.PipelineControl import *


@files([Paths.genomespace_refseq_counts(sample)
        for sample in Options.ALLCLIP_SAMPLES],
       Paths.clip_enrichment_summary)
@follows(DerivedDatabaseBuilding.quantitate_refseq_in_gspace)
def clip_refseq_enrichment_statistics(inputfiles, outputfile):
    inputfilesarg = ' '.join('%s:%s' % (
        smp, Paths.genomespace_refseq_counts(smp))
        for smp in Options.ALLCLIP_SAMPLES)
    refsample = Options.CLIPCTL_SAMPLES[0] # XXX fix this to support multiple controls
    clipsamples = ','.join(Options.CLIP_SAMPLES)

    runproc("""
        $STATS_CLIP_NRREFSEQ_ENRICHED \
            $nr_refseq_db $refsample $clipsamples $inputfilesarg \
            > $outputfile""", outputfile)


@files(for_each([Paths.fulltag_annotation, Paths.fulltag_prealn_sam],
                Paths.total_read_class_stats, Options.ALL_SAMPLES))
@follows(SequenceAnnotation.summarize_annotations)
@follows(ContaminantFilter.align_to_known_contaminants)
def total_read_classification_statistics(inputfiles, outputfile, sample):
    annofile, prealnfile = inputfiles

    with TemporaryFile() as tmpstats:
        runproc("$STATS_READ_CLASS_PROPORTION $annofile > $tmpstats")
        runproc("$STATS_READ_CLASS_ADD_FILTERED $tmpstats $prealnfile "
                "> $outputfile", outputfile)


@files(for_each(Paths.fulltag_collapsed_reads,
                Paths.read_length_stats, Options.ALL_SAMPLES))
@follows(SequenceProcessing.fulltag_collapse)
def read_length_statistics(inputfile, outputfile, sample):
    runproc("""
        $STATS_LENGTH_DIST $inputfile > $outputfile""", outputfile)


def tasks():
    return [
        clip_refseq_enrichment_statistics,
        total_read_classification_statistics,
        read_length_statistics,
    ]
