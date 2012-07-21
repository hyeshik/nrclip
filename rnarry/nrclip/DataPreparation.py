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
import tempfile
import shutil
from rnarry.nrclip import Paths
from rnarry.nrclip.PipelineControl import *

@files(Paths.early_filter_fasta, Paths.early_filter_index_check)
def generate_gsnap_early_filter_index(inputfile, outcheck):
    command = '%s -d %s -D %s %s' % (
        Paths.GMAP_BUILD_CMD, Paths.early_filter_prefix,
        Paths.SETTINGS_DIR, inputfile)
    runproc(command)

GENOME_SEQ_URL = 'http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz'
@files(None, Paths.genome_fasta_zipped)
def download_zipped_genome_sequence(inputfile, outputfile):
    command = '%s -O %s %s' % (Paths.WGET_CMD, outputfile, GENOME_SEQ_URL)
    runproc(command)

@files(Paths.genome_fasta_zipped, Paths.genome_fasta)
@follows(download_zipped_genome_sequence)
def extract_genome_sequence(inputfile, outputfile):
    tdir = tempfile.mkdtemp(dir=Paths.TMP_DIR)
    try:
        runproc('%s -C %s -xzf %s' % (Paths.TAR_CMD, tdir, inputfile))
        runproc('cd %s && cat *.fa > %s' % (tdir, outputfile))
    finally:
        shutil.rmtree(tdir)

def tasks():
    return [
        generate_gsnap_early_filter_index,
        download_zipped_genome_sequence,
        extract_genome_sequence,
    ]

