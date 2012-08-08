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


@files(for_each([Paths.fulltag_masked_reads, Paths.reftranscriptome_gmap_index],
                Paths.fulltag_transcriptome_alignment_gmap,
                Paths.FULLTAG_SAMPLES,
                [Options.FULLTAG_TRANSCRIPTOME_MISMATCHES]))
@follows(DerivedDatabaseBuilding.generate_gsnap_transcriptome_index)
@follows(SequenceMasking.produce_masked_fasta)
@jobs_limit(1, 'exclusive')
def align_to_transcriptome_gmap(inputfiles, outputfile, sample, mismatches):
    reads, indexcheck = inputfiles

    runproc("""
        $GSNAP -D $REFTRANSCRIPTOME_DIR -d $reftranscriptome_dbname -O -B 4 \
            --terminal-threshold=999999 --trim-mismatch-score=0 \
            -m $mismatches -t $NUM_THREADS $reads | \
        $GZIP -c - > $outputfile""", outputfile)


@files(for_each(Paths.fulltag_transcriptome_alignment_sam,
                Paths.fulltag_transcriptomic_besthits_sam,
                Paths.FULLTAG_SAMPLES,
                [Options.FULLTAG_POSTPROC_ALLOWED_MISMATCHES]))
@follows(align_to_transcriptome_sam)
def resolve_transcriptomic_multihits_sam(inputfile, outputfile, sample, mismatches):
    runproc("""
        $ZCAT $inputfile | \
        $AWK -F'\t' '
            /^@/ { print $$0; }
            /^[^@]/ {
                if (and($$2, 16) == 0)
                    print $$0;
            }' | \
        $SAM_MULTIHIT_RESOLVE $mismatches | \
        $GZIP -c - > $outputfile""", outputfile)


@files(for_each(Paths.fulltag_transcriptome_alignment_gmap,
                Paths.fulltag_transcriptomic_besthits_gmap,
                Paths.FULLTAG_SAMPLES,
                [Options.FULLTAG_POSTPROC_ALLOWED_MISMATCHES]))
@follows(align_to_transcriptome_gmap)
def resolve_transcriptomic_multihits_gmap(inputfile, outputfile, gmapple, mismatches):
    runproc("""
        $ZCAT $inputfile | \
        $GMAP_MULTIHIT_RESOLVE $mismatches -tr | \
        $GZIP -c - > $outputfile""", outputfile)


@files(for_each(Paths.fulltag_transcriptomic_besthits_sam,
                Paths.tspace_read_database,
                Paths.FULLTAG_SAMPLES))
@follows(resolve_transcriptomic_multihits_sam)
def build_tspace_read_database(inputfile, outputfile, sample):
    runproc("""
        $PREPARE_FLATDATA_FROM_SAM $inputfile -tr | \
        $SORT -k3,4 -k5,5n | \
        $BUILD_POSITIONALDB_TRANSCRIPTOME $nr_refseq_db $outputfile""",
        outputfile)


def tasks():
    return [
        align_to_transcriptome_sam,
        align_to_transcriptome_gmap,
        resolve_transcriptomic_multihits_sam,
        resolve_transcriptomic_multihits_gmap,
        build_tspace_read_database,
    ]

