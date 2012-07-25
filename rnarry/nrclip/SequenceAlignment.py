#
# rnarry.nrclip.SequenceAlignment
#  - aligns sequences to genome or transcriptome
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
from rnarry.nrclip import Paths, Options, ContaminantFilter, DataPreparation
from rnarry.nrclip.PipelineControl import *


@files(for_each_sample(Paths.fulltag_filtered_reads,
                       Paths.fulltag_genome_alignment_sam,
                       Paths.ALL_SAMPLES))
@follows(DataPreparation.generate_gsnap_genome_index)
@follows(DataPreparation.build_splicesites_index)
@follows(ContaminantFilter.fulltag_filter_contaminant)
@jobs_limit(1, 'exclusive')
def fulltag_genome_alignment_sam(inputfile, outputfile, sample):
    runproc("""
        $GSNAP -D $EXTERNAL_DIR -d $genome_prefix -O -B 4 -A sam \
            --terminal-threshold=9999 -s $splice_index \
            -m $FULLTAG_GENOME_MISMATCHES -t $NUM_THREADS $inputfile | \
        $GZIP -c - > $outputfile""", outputfile)


@files(for_each_sample(Paths.fulltag_genome_alignment_sam,
                       Paths.fulltag_genome_besthits,
                       Paths.ALL_SAMPLES))
@follows(fulltag_genome_alignment_sam)
def fulltag_genome_resolve_multihit(inputfile, outputfile, sample):
    runproc("""
        $ZCAT $inputfile | \
        $SAM_MULTIHIT_RESOLVE $GENOMEALN_POSTPROC_ALLOWED_MISMATCHES | \
        $GZIP -c - > $outputfile""", outputfile)


@files(for_each_sample(Paths.fulltag_genome_alignment_sam,
                       Paths.fulltag_genome_alignment_unsorted_bam,
                       Paths.ALL_SAMPLES))
@follows(fulltag_genome_alignment_sam)
def fulltag_convert_primary_alignment_to_bam(inputfile, outputfile, sample):
    runproc("$SAMTOOLS view -bS -o $outputfile $inputfile")


@files(for_each_sample(Paths.shorttag_filtered_reads,
                       Paths.shorttag_genome_alignment_sam,
                       Paths.SHORTTAG_SAMPLES))
@follows(DataPreparation.generate_gsnap_genome_index)
@follows(DataPreparation.build_splicesites_index)
@follows(ContaminantFilter.shorttag_filter_contaminant)
@jobs_limit(1, 'exclusive')
def shorttag_genome_alignment_sam(inputfile, outputfile, sample):
    runproc("""
        $GSNAP -D $EXTERNAL_DIR -d $genome_prefix -O -B 4 -A sam \
            --terminal-threshold=9999 -s $splice_index \
            -m $FULLTAG_GENOME_MISMATCHES -t $NUM_THREADS $inputfile | \
        $GZIP -c - > $outputfile""", outputfile)


@files(for_each_sample(Paths.shorttag_genome_alignment_sam,
                       Paths.shorttag_genome_alignment_unsorted_bam,
                       Paths.SHORTTAG_SAMPLES))
@follows(shorttag_genome_alignment_sam)
def shorttag_convert_primary_alignment_to_bam(inputfile, outputfile, sample):
    runproc("$SAMTOOLS view -bS -o $outputfile $inputfile")


def tasks():
    return [
        fulltag_genome_alignment_sam,
        fulltag_genome_resolve_multihit,
        fulltag_convert_primary_alignment_to_bam,
        shorttag_genome_alignment_sam,
        shorttag_convert_primary_alignment_to_bam,
    ]

