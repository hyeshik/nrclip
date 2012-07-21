#
# rnarry.nrclip.Paths
#  - Sets paths to data files
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

import os

TOPDIR = os.getcwd()

def rel(path):
    return os.path.join(TOPDIR, path)


# =================================
# Configurations for external tools

GSNAP_CMD = 'gsnap'
GMAP_BUILD_CMD = 'gmap_build'
WGET_CMD = 'wget'
TAR_CMD = 'tar'


# ==============
# Subdirectories

EXTERNAL_DIR = 'external'
TMP_DIR = 'tmp'


# ===============================
# Settings and external resources

early_filter_prefix = 'early-filter'
early_filter_fasta = rel('%s/%s.fa' % (EXTERNAL_DIR, early_filter_prefix))
early_filter_index_check = rel('%s/%s/%s.ref153positions' % (
                    EXTERNAL_DIR, early_filter_prefix, early_filter_prefix))
genome_prefix = 'mm9'
genome_fasta = rel('%s/%s.fa' % (EXTERNAL_DIR, genome_prefix))
genome_index_check = rel('%s/%s/%s.ref12123positions' % (
                         EXTERNAL_DIR, genome_prefix, genome_prefix))
genome_fasta_zipped = rel('%s/%s.tar.gz' % (EXTERNAL_DIR, genome_prefix))


