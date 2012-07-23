#
# rnarry.nrclip.DataPreparation
#  - Downloads and preprocesses external data and settings.
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
from rnarry.nrclip import Paths
from rnarry.nrclip.PipelineControl import *


@files(Paths.early_filter_fasta, Paths.early_filter_index_check)
def generate_gsnap_early_filter_index(inputfile, outcheck):
    runproc('$GMAP_BUILD -d $early_filter_prefix '
            '-D $EXTERNAL_DIR $inputfile', outcheck)


@files(None, Paths.genome_fasta_zipped)
def download_zipped_genome_sequence(inputfile, outputfile):
    runproc('$WGET -O $outputfile $GENOME_SEQ_URL', outputfile)


@files(Paths.genome_fasta_zipped, Paths.genome_fasta)
@follows(download_zipped_genome_sequence)
def extract_genome_sequence(inputfile, outputfile):
    with TemporaryDirectory() as tdir:
        runproc('$TAR -C $tdir -xzf $inputfile')
        runproc('cd $tdir && cat *.fa > $outputfile', outputfile)


@files(Paths.genome_fasta, Paths.genome_index_check)
@follows(extract_genome_sequence)
def generate_gsnap_genome_index(inputfile, outcheck):
    runproc('$GMAP_BUILD -d $genome_prefix -D $EXTERNAL_DIR -k 12 '
            '$inputfile', outcheck)


@files(None, Paths.refgene_ucsc)
def download_refgene(inputfile, outputfile):
    runproc('$WGET -O $outputfile $REFGENE_URL', outputfile)


@files(None, Paths.knowngene_ucsc)
def download_knowngene(inputfile, outputfile):
    runproc('$WGET -O $outputfile $KNOWNGENE_URL', outputfile)


@files([Paths.refgene_ucsc, Paths.knowngene_ucsc], Paths.splice_index)
@follows(download_refgene)
@follows(download_knowngene)
def build_splicesites_index(inputfiles, outputfile):
    refgene, knowngene = inputfiles

    runproc("""
        ($ZCAT $refgene | $PSL_SPLICESITES -s 1; \
         $ZCAT $knowngene | $PSL_SPLICESITES) | \
        $IIT_STORE -o $outputfile""", outputfile)


def tasks():
    return [
        generate_gsnap_early_filter_index,
        download_zipped_genome_sequence,
        extract_genome_sequence,
        generate_gsnap_genome_index,
        download_refgene,
        download_knowngene,
        build_splicesites_index,
    ]
