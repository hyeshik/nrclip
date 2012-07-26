#
# rnarry.nrclip.DerivedDatabaseBuilding
#  - Builds databases derived from analyses of experimental data
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
from rnarry.nrclip import Paths, Options, SequenceMasking, DataPreparation
from rnarry.nrclip.PipelineControl import *


@files(for_each(Paths.fulltag_genome_besthits,
                Paths.genomespace_read_database,
                Paths.ALL_SAMPLES))
@follows(SequenceMasking.resolve_genomic_multihits)
def build_genomespace_read_database(inputfile, outputfile, sample):
    outputdir = os.path.dirname(outputfile)
    runproc("""
        $PREPARE_FLATDATA_FROM_SAM $inputfile | \
        $SORT -k1,1 -k3,4 -k5,5n | \
        $BUILD_POSITIONALDB_GENOME $outputdir""", outputfile)


@files(for_each(Paths.genomespace_read_database,
                Paths.genomespace_refseq_counts,
                Paths.ALL_SAMPLES))
@follows(build_genomespace_read_database)
@follows(DataPreparation.build_nonredundant_refseq_database)
def quantitate_refseq_in_gspace(inputfile, outputfile, sample):
    gspacedir = os.path.dirname(inputfile)
    runproc(
        '$COUNT_REFSEQ_IN_GSPACE $outputfile $nr_refseq_db $gspacedir',
        outputfile)


def tasks():
    return [
        build_genomespace_read_database,
        quantitate_refseq_in_gspace,
    ]
