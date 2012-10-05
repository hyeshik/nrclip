#
# rnarry.nrclip.RPFStatistics
#  - Produces ribosome profiling-specific statistics
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
    Paths, Options, DataPreparation, TranscriptomeAnalysis)
from rnarry.nrclip.PipelineControl import *


@files(for_each([Paths.nr_refseq_db, Paths.tspace_read_database],
                Paths.cds_read_count_table,
                Options.SHORTTAG_SAMPLES))
@follows(DataPreparation.build_nonredundant_refseq_database)
@follows(TranscriptomeAnalysis.build_tspace_read_database)
def count_cds_reads(inputfiles, outputfile, sample):
    nr_refseq_db, tspace = inputfiles

    runproc("""
        $TSPACE_COUNT_CDS $nr_refseq_db $tspace > $outputfile""", outputfile)


def tasks():
    return [
        count_cds_reads,
    ]
