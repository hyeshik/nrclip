#
# rnarry.nrclip.ErrorAnalysis
#  - Error profile analysis for CLIP-seq samples
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
from rnarry.nrclip import Paths, Options, TranscriptomeAnalysis
from rnarry.nrclip.PipelineControl import *
from rnarry.sequtils import get_first_sequence_length


@files(for_each([Paths.fulltag_transcriptomic_besthits_gmap,
                 Paths.original_sequence_reads],
                Paths.error_profile_read_level,
                Paths.ALLCLIP_SAMPLES))
@follows(TranscriptomeAnalysis.resolve_transcriptomic_multihits_gmap)
def count_alignment_errors(inputfiles, outputfile, sample):
    mapping, origseqs = inputfiles
    readlength = get_first_sequence_length(origseqs)

    runproc("""
        $GMAP_ERROR_PROFILE $mapping $outputfile $readlength""", outputfile)


def tasks():
    return [
        count_alignment_errors,
    ]

