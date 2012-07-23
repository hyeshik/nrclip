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
from rnarry.nrclip.Options import *


TOPDIR = os.getcwd()

def rel(path):
    return os.path.join(TOPDIR, path)

def relfmt(path):
    return lambda *v: os.path.join(TOPDIR, path % tuple(v))


# =====================
# External tool aliases

#- Generic UNIX command line tools
AWK = 'awk'
CUT = 'cut'
GREP = 'grep'
GZIP = 'pigz'
GZIP_LT = 'gzip' # for low throughput outputs
SORT = 'sort'
TAR = 'tar'
UNIQ = 'uniq'
WGET = 'wget'
ZCAT = 'zcat'

#- A. Gordon's FASTX_Toolkit
FASTX_CLIPPER = 'fastx_clipper'
FASTQ_QUALITY_TRIMMER = 'fastq_quality_trimmer'
FASTQ_QUALITY_FILTER = 'fastq_quality_filter'
FASTX_COLLAPSER = 'fastx_collapser'
FASTX_TRIMMER = 'fastx_trimmer'
FASTX_ARTIFACTS_FILTER = 'fastx_artifacts_filter'

#- Bioinformatic tools
GSNAP = 'gsnap' # from gmap/gsnap
GMAP_BUILD = 'gmap_build' # from gmap/gsnap
FASOMERECORDS = 'faSomeRecords' # from Jim Kent's ucscgb
PSL_SPLICESITES = 'psl_splicesites' # from Jim Kent's ucscgb
IIT_STORE = 'iit_store' # from gmap/gsnap


# =================================
# External resource URLs

GENOME_SEQ_URL = 'http://hgdownload.cse.ucsc.edu/goldenPath/%s/bigZips/chromFa.tar.gz' % GENOME
REFGENE_URL = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/%s/database/refGene.txt.gz' % GENOME
KNOWNGENE_URL = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/%s/database/knownGene.txt.gz' % GENOME


# ==============
# Subdirectories

EXTERNAL_DIR = 'external'
ORIGREAD_DIR = 'sequences'
GENOMEALN_DIR = 'genomealn'
PREALN_DIR = 'prealns'
TMP_DIR = 'tmp'


# ===============================
# Settings and external resources

early_filter_prefix = 'early-filter'
early_filter_fasta = rel('%s/%s.fa' % (EXTERNAL_DIR, early_filter_prefix))
early_filter_index_check = rel('%s/%s/%s.ref153positions' % (
                    EXTERNAL_DIR, early_filter_prefix, early_filter_prefix))
genome_prefix = GENOME
genome_fasta = rel('%s/%s.fa' % (EXTERNAL_DIR, genome_prefix))
genome_index_check = rel('%s/%s/%s.ref123positions' % (
                         EXTERNAL_DIR, genome_prefix, genome_prefix))
genome_fasta_zipped = rel('%s/%s.tar.gz' % (EXTERNAL_DIR, genome_prefix))
refgene_ucsc = rel(EXTERNAL_DIR + '/refGene.txt.gz')
knowngene_ucsc = rel(EXTERNAL_DIR + '/knownGene.txt.gz')
splice_index = rel(EXTERNAL_DIR + '/%s/%s.maps/%s.splicesites.iit' % (
                    genome_prefix, genome_prefix, genome_prefix))


# =========================
# Basic sequence processing

original_sequence_reads = relfmt(ORIGREAD_DIR + '/%s.fq.gz')
fulltag_quality_filtered_reads = relfmt(ORIGREAD_DIR + '/%s-hq.fq.gz')
fulltag_collapsed_reads = relfmt(ORIGREAD_DIR + '/%s-total.fa')
shorttag_trimmed_reads = relfmt(ORIGREAD_DIR + '/%s-trimmed.fq.gz')
shorttag_tags = relfmt(ORIGREAD_DIR + '/%s-tags.fa')


# ==================
# Contaminant filter

fulltag_prealn_sam = relfmt(PREALN_DIR + '/full-%s.sam.gz')
shorttag_prealn_sam = relfmt(PREALN_DIR + '/short-%s.sam.gz')

fulltag_filtered_reads = relfmt(PREALN_DIR + '/full-%s-clean.fa')
shorttag_filtered_reads = relfmt(PREALN_DIR + '/short-%s-clean.fa')


# ===================
# Sequence alignments

fulltag_genome_alignment_sam = relfmt(GENOMEALN_DIR + '/%s.sam.gz')

