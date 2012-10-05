#
# rnarry.nrclip.SequenceAnnotation
#  - Sequence annotation, classifcation, and related statistics
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
from rnarry.nrclip import Paths, Options, DataPreparation, SequenceAlignment
from rnarry.nrclip.PipelineControl import *


@files(for_each([Paths.fulltag_genome_alignment_unsorted_bam,
                 Paths.compiled_catalog],
                Paths.fulltag_primary_annotation,
                Paths.ALL_SAMPLES) +
       for_each([Paths.shorttag_genome_alignment_unsorted_bam,
                 Paths.compiled_catalog],
                Paths.shorttag_primary_annotation,
                Paths.SHORTTAG_SAMPLES))
@follows(DataPreparation.compile_catalogs)
@follows(SequenceAlignment.convert_primary_alignment_to_bam)
def annotate_sequences(inputfiles, outputfile, sample):
    inputfile, catalog = inputfiles
    runproc("""
        $INTERSECTBED -bed -abam $inputfile -b $catalog -wa -wb | \
        $GZIP -c - > $outputfile""", outputfile)


@files(for_each([Paths.fulltag_genome_alignment_sam,
                 Paths.fulltag_primary_annotation],
                Paths.fulltag_annotation,
                Paths.ALL_SAMPLES) +
       for_each([Paths.shorttag_genome_alignment_sam,
                 Paths.shorttag_primary_annotation],
                Paths.shorttag_annotation,
                Paths.SHORTTAG_SAMPLES))
@follows(annotate_sequences)
def summarize_annotations(inputfiles, outputfile, sample):
    saminput, bedinput = inputfiles
    runproc("""
        $SUMMARIZE_ANNOTATIONS $saminput $bedinput | $GZIP -c - > $outputfile
        """, outputfile)


#@files(for_each(Paths.shorttag_annotation,
#                Paths.shorttag_mRNA_reads_count,
#                Paths.SHORTTAG_SAMPLES))
#@follows(summarize_annotations)
#def count_mRNA_reads(inputfile, outputfile, sample):
#    runproc("""
#        $ZCAT $inputfile |
#        $AWK -F'\t' 'BEGIN {cnt=0}
#            { if ($$3 == "CDS" || $$3 == "UTR5" || $$3 == "UTR3") {
#                split($$1, tok, "-");
#                cnt += tok[2];
#              }
#            } END { print cnt }' > $outputfile
#        """, outputfile)


def tasks():
    return [
        annotate_sequences,
        summarize_annotations,
#        count_mRNA_reads,
    ]

