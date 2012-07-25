#
# rnarry.nrclip.SequenceMasking
#  - Masks rRNA and tRNA sequences from reads or alignments
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
from rnarry.nrclip import Paths, Options, SequenceAnnotation
from rnarry.nrclip.PipelineControl import *


@files(for_each(Paths.fulltag_primary_annotation,
                Paths.fulltag_masked_readids,
                Paths.ALL_SAMPLES) +
       for_each(Paths.shorttag_primary_annotation,
                Paths.shorttag_masked_readids,
                Paths.SHORTTAG_SAMPLES))
@follows(SequenceAnnotation.annotate_sequences)
def make_list_of_masked_sequences(inputfile, outputfile, sample):
    runproc("$ZGREP '	[rt]RNA' $inputfile | $CUT -f4 | $UNIQ > $outputfile",
            outputfile)


@files(for_each([Paths.fulltag_genome_alignment_sam,
                 Paths.fulltag_masked_readids],
                Paths.fulltag_masked_genome_alignments,
                Paths.ALL_SAMPLES) +
       for_each([Paths.shorttag_genome_alignment_sam,
                 Paths.shorttag_masked_readids],
                Paths.shorttag_masked_genome_alignments,
                Paths.SHORTTAG_SAMPLES))
@follows(make_list_of_masked_sequences)
def produce_masked_sam(inputfile, outputfile, sample):
    alnfile, maskinglist = inputfile
    runproc("$SAM_ID_FILTER $alnfile $maskinglist | $GZIP -c - > $outputfile",
            outputfile)


@files(for_each([Paths.fulltag_filtered_reads, Paths.fulltag_masked_readids],
                Paths.fulltag_masked_reads,
                Paths.ALL_SAMPLES) +
       for_each([Paths.shorttag_filtered_reads, Paths.shorttag_masked_readids],
                Paths.shorttag_masked_reads,
                Paths.SHORTTAG_SAMPLES))
@follows(make_list_of_masked_sequences)
def produce_masked_fasta(inputfile, outputfile, sample):
    seqfile, maskinglist = inputfile
    runproc("$FASOMERECORDS -exclude $seqfile $maskinglist $outputfile",
            outputfile)



@files(for_each(Paths.fulltag_masked_genome_alignments,
                Paths.fulltag_genome_besthits,
                Paths.ALL_SAMPLES,
                [Options.FULLTAG_POSTPROC_ALLOWED_MISMATCHES]) +
       for_each(Paths.shorttag_masked_genome_alignments,
                Paths.shorttag_genome_besthits,
                Paths.SHORTTAG_SAMPLES,
                [Options.SHORTTAG_POSTPROC_ALLOWED_MISMATCHES]))
@follows(produce_masked_sam)
def resolve_genomic_multihits(inputfile, outputfile, sample, mismatches):
    runproc("""
        $ZCAT $inputfile | $SAM_MULTIHIT_RESOLVE $mismatches | \
        $GZIP -c - > $outputfile""", outputfile)


def tasks():
    return [
        make_list_of_masked_sequences,
        produce_masked_sam,
        produce_masked_fasta,
        resolve_genomic_multihits,
    ]

