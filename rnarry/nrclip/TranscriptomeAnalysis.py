#
# rnarry.nrclip.TranscriptomeAnalysis
#  - alignment and analyses against transcriptome
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
from rnarry.nrclip import (
    Paths, Options, DerivedDatabaseBuilding, SequenceMasking)
from rnarry.nrclip.PipelineControl import *


@files(for_each([Paths.fulltag_masked_reads, Paths.reftranscriptome_gmap_index],
                Paths.fulltag_transcriptome_alignment_sam,
                Paths.FULLTAG_SAMPLES,
                [Options.FULLTAG_TRANSCRIPTOME_MISMATCHES]))
@follows(DerivedDatabaseBuilding.generate_gsnap_transcriptome_index)
@follows(SequenceMasking.produce_masked_fasta)
@jobs_limit(1, 'exclusive')
def align_to_transcriptome_sam(inputfiles, outputfile, sample, mismatches):
    reads, indexcheck = inputfiles

    runproc("""
        $GSNAP -D $REFTRANSCRIPTOME_DIR -d $reftranscriptome_dbname -O -B 4 \
            -A sam --terminal-threshold=9999 --trim-mismatch-score=0 \
            -m $mismatches -t $NUM_THREADS $reads | \
        $GZIP -c - > $outputfile""", outputfile)


#@files(for_each(Paths.fulltag_genome_alignment_sam,
#                Paths.fulltag_genome_alignment_unsorted_bam,
#                Paths.ALL_SAMPLES) +
#       for_each(Paths.shorttag_genome_alignment_sam,
#                Paths.shorttag_genome_alignment_unsorted_bam,
#                Paths.SHORTTAG_SAMPLES))
#@follows(align_to_genome_sam)
#def convert_primary_alignment_to_bam(inputfile, outputfile, sample):
#    runproc("$SAMTOOLS view -bS -o $outputfile $inputfile")


def tasks():
    return [
        align_to_transcriptome_sam,
        #convert_primary_alignment_to_bam,
    ]

