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
GZIP = 'gzip'
SED = 'sed'
SORT = 'sort'
TAR = 'tar'
UNIQ = 'uniq'
WGET = 'wget'
ZCAT = 'zcat'
ZGREP = 'zgrep'

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
IIT_STORE = 'iit_store' # from gmap/gsnap
INTERSECTBED = 'intersectBed' # from bedtools
SAMTOOLS = 'samtools' # from samtools

#- Jim Kent's UCSC Genome Browser tools
FASOMERECORDS = 'faSomeRecords'
PSL_SPLICESITES = 'psl_splicesites'
FA2TWOBIT = 'faToTwoBit'
BLAT = 'blat'
PSL_CDNAFILTER = 'pslCDnaFilter'
PSL2BED = 'pslToBed'

#- Local scripts
SAM_MULTIHIT_RESOLVE = 'python scripts/sam-multihit-resolve.py'
BUILD_REFSEQ_INDEX = 'python scripts/build-refseq-index.py'
BUILD_RFAM_INDEX = 'sh scripts/build-rfam-index.sh'


# =================================
# External resource URLs

GENOME_SEQ_URL = 'http://hgdownload.cse.ucsc.edu/goldenPath/%s/bigZips/chromFa.tar.gz' % GENOME
REFGENE_URL = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/%s/database/refGene.txt.gz' % GENOME
KNOWNGENE_URL = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/%s/database/knownGene.txt.gz' % GENOME
MIRBASE_URL = 'ftp://mirbase.org/pub/mirbase/18/genomes/%s.gff2' % SPECIES3
REPEATMASKER_URL = 'http://hgdownload.cse.ucsc.edu/goldenPath/%s/database/chr%s_rmsk.txt.gz'
RFAM_FULL_URL = 'ftp://ftp.sanger.ac.uk/pub/databases/Rfam/CURRENT/Rfam.full.gz'
RFAM_SEQUENCE_URL = 'ftp://ftp.sanger.ac.uk/pub/databases/Rfam/CURRENT/Rfam.fasta.gz'
GTRNADB_URL = 'http://hgdownload.cse.ucsc.edu/goldenPath/%s/database/tRNAs.txt.gz' % GENOME


# ==============
# Subdirectories

EXTERNAL_DIR = 'external'
ORIGREAD_DIR = 'sequences'
GENOMEALN_DIR = 'genomealn'
ANNOTATIONS_DIR = 'annotations'
DOWNLOAD_DIR = 'external/downloaded'
PREALN_DIR = 'prealns'
TMP_DIR = 'tmp'

ALL_SUBDIRS = [
    EXTERNAL_DIR, ORIGREAD_DIR, GENOMEALN_DIR, DOWNLOAD_DIR,
    PREALN_DIR, TMP_DIR, ANNOTATIONS_DIR
]


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
genome_twobit = rel('%s/%s.2bit' % (EXTERNAL_DIR, genome_prefix))
refgene_ucsc = rel(EXTERNAL_DIR + '/refGene.txt.gz')
knowngene_ucsc = rel(EXTERNAL_DIR + '/knownGene.txt.gz')
splice_index = rel(EXTERNAL_DIR + '/%s/%s.maps/%s.splicesites.iit' % (
                    genome_prefix, genome_prefix, genome_prefix))

mirbase_catalog = rel(EXTERNAL_DIR + '/mirbase.bed.gz')
refseq_catalog = rel(EXTERNAL_DIR + '/refseq.bed.gz')
repeatmasker_original = relfmt(DOWNLOAD_DIR + '/rmsk-chr%s.txt.gz')
repeatmasker_catalog = rel(EXTERNAL_DIR + '/repeatmasker.bed.gz')
rfam_fasta = rel(DOWNLOAD_DIR + '/Rfam.fasta.gz')
rfam_original = rel(DOWNLOAD_DIR + '/Rfam.full.gz')
rfam_catalog = rel(EXTERNAL_DIR + '/rfam.bed.gz')
trna_catalog = rel(EXTERNAL_DIR + '/trna.bed.gz')
compiled_catalog = rel(EXTERNAL_DIR + '/all-annotations.bed.gz')


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
fulltag_genome_alignment_unsorted_bam = relfmt(GENOMEALN_DIR + '/%s.unsorted.bam')
fulltag_genome_besthits = relfmt(GENOMEALN_DIR + '/%s-besthits.sam.gz')


# ===========
# Annotations

fulltag_primary_annotation = relfmt(ANNOTATIONS_DIR + '/%s.anno.gz')

