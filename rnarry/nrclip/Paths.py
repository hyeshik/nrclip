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
AWK = 'gawk'
CUT = 'cut'
ENV = 'env'
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
FASIZE = 'faSize'

#- In-house scripts
PYTHON = 'python'
SAM_MULTIHIT_RESOLVE = PYTHON + ' scripts/sam-multihit-resolve.py'
GMAP_MULTIHIT_RESOLVE = PYTHON + ' scripts/gmap-multihit-resolve.py'
BUILD_REFSEQ_INDEX = PYTHON + ' scripts/build-refseq-index.py'
BUILD_RFAM_INDEX = 'sh scripts/build-rfam-index.sh'
SUMMARIZE_ANNOTATIONS = PYTHON + ' scripts/summarize-annotations.py'
SAM_ID_FILTER = PYTHON + ' scripts/sam-id-filter.py'
PREPARE_FLATDATA_FROM_SAM = PYTHON + ' scripts/prepare-flatdata-from-sam.py'
BUILD_POSITIONALDB_GENOME = PYTHON + ' scripts/build-positionaldb-genome.py'
BUILD_POSITIONALDB_TRANSCRIPTOME = PYTHON + ' scripts/build-positionaldb-transcriptome.py'
BUILD_NONREDUNDANT_REFSEQ = PYTHON + ' scripts/build-nonredundant-refseq.py'
NRREFSEQ2BED = PYTHON + ' scripts/nrrefseq2bed.py'
COUNT_REFSEQ_IN_GSPACE = PYTHON + ' scripts/count-refseq-in-gspace.py'
REFSEQCNT_PICK_EXPRESSED = PYTHON + ' scripts/refseqcnt-pick-expressed.py'
GENFASTA_MUTATED_TRANSCRIPTOME = PYTHON + ' scripts/genfasta-mutated-transcriptome.py'
STATS_CLIP_NRREFSEQ_ENRICHED = PYTHON + ' scripts/clip-nrrefseq-enriched.py'
STATS_READ_CLASS_PROPORTION = PYTHON + ' scripts/stats-read-class-proportion.py'
STATS_READ_CLASS_ADD_FILTERED = 'sh scripts/stats-classstat-addprealns.sh'
GMAP_ERROR_PROFILE = PYTHON + ' scripts/gmap-error-profile.py'


# =================================
# External resource URLs

GENOME_SEQ_URL = 'http://hgdownload.cse.ucsc.edu/goldenPath/%s/bigZips/chromFa.tar.gz' % GENOME
REFGENE_URL = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/%s/database/refGene.txt.gz' % GENOME
REFFLAT_URL = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/%s/database/refFlat.txt.gz' % GENOME
REFLINK_URL = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/%s/database/refLink.txt.gz' % GENOME
KNOWNGENE_URL = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/%s/database/knownGene.txt.gz' % GENOME
MIRBASE_URL = 'ftp://mirbase.org/pub/mirbase/18/genomes/%s.gff2' % SPECIES3
REPEATMASKER_URL = 'http://hgdownload.cse.ucsc.edu/goldenPath/%s/database/chr%s_rmsk.txt.gz'
RFAM_FULL_URL = 'ftp://ftp.sanger.ac.uk/pub/databases/Rfam/CURRENT/Rfam.full.gz'
RFAM_SEQUENCE_URL = 'ftp://ftp.sanger.ac.uk/pub/databases/Rfam/CURRENT/Rfam.fasta.gz'
GTRNADB_URL = 'http://hgdownload.cse.ucsc.edu/goldenPath/%s/database/tRNAs.txt.gz' % GENOME


# ==============
# Subdirectories

EXTERNAL_DIR = 'external'
SEQUENCES_DIR = 'sequences'
GENOMEALN_DIR = 'genomealn'
ANNOTATIONS_DIR = 'annotations'
DOWNLOAD_DIR = 'external/downloaded'
DERIVEDBASES_DIR = 'derived'
REFTRANSCRIPTOME_DIR = 'derived/reference'
STATISTICS_DIR = 'stats'
TRANSCRIPTOMEALN_DIR = 'transcriptomealn'
PREALN_DIR = 'prealns'
ERRORANALYSIS_DIR = 'erroranalysis'
TMP_DIR = 'tmp'

ALL_SUBDIRS = [
    EXTERNAL_DIR, SEQUENCES_DIR, GENOMEALN_DIR, DOWNLOAD_DIR,
    PREALN_DIR, TMP_DIR, ANNOTATIONS_DIR, DERIVEDBASES_DIR,
    REFTRANSCRIPTOME_DIR, TRANSCRIPTOMEALN_DIR, STATISTICS_DIR,
    ERRORANALYSIS_DIR,
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
refflat_ucsc = rel(EXTERNAL_DIR + '/refFlat.txt.gz')
reflink_ucsc = rel(EXTERNAL_DIR + '/refLink.txt.gz')
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

nr_refseq_db = rel(EXTERNAL_DIR + '/nrRefSeq.db')
nr_refseq_list = rel(EXTERNAL_DIR + '/nrRefSeq.list')
nr_refseq_genome_bed = rel(EXTERNAL_DIR + '/nrRefSeq-genome.bed.gz')


# =========================
# Basic sequence processing

original_sequence_reads = relfmt(SEQUENCES_DIR + '/%s.fq.gz')
fulltag_quality_filtered_reads = relfmt(SEQUENCES_DIR + '/%s-hq.fq.gz')
fulltag_collapsed_reads = relfmt(SEQUENCES_DIR + '/%s-total.fa')
shorttag_trimmed_reads = relfmt(SEQUENCES_DIR + '/%s-trimmed.fq.gz')
shorttag_tags = relfmt(SEQUENCES_DIR + '/%s-tags.fa')


# ==================
# Contaminant filter

fulltag_prealn_sam = relfmt(PREALN_DIR + '/full-%s.sam.gz')
shorttag_prealn_sam = relfmt(PREALN_DIR + '/short-%s.sam.gz')

fulltag_filtered_reads = relfmt(SEQUENCES_DIR + '/full-%s-filtered.fa')
shorttag_filtered_reads = relfmt(SEQUENCES_DIR + '/short-%s-filtered.fa')


# ===================
# Sequence alignments

fulltag_genome_alignment_sam = relfmt(GENOMEALN_DIR + '/%s.sam.gz')
fulltag_genome_alignment_unsorted_bam = relfmt(GENOMEALN_DIR + '/%s.unsorted.bam')
fulltag_genome_besthits = relfmt(GENOMEALN_DIR + '/%s-besthits.sam.gz')
fulltag_masked_genome_alignments = relfmt(GENOMEALN_DIR + '/%s-masked.sam.gz')

shorttag_genome_alignment_sam = relfmt(GENOMEALN_DIR + '/short-%s.sam.gz')
shorttag_genome_alignment_unsorted_bam = relfmt(GENOMEALN_DIR + '/short-%s.unsorted.bam')
shorttag_genome_besthits = relfmt(GENOMEALN_DIR + '/short-%s-besthits.sam.gz')
shorttag_masked_genome_alignments = relfmt(GENOMEALN_DIR + '/short-%s-masked.sam.gz')

# ===========
# Annotations

fulltag_primary_annotation = relfmt(ANNOTATIONS_DIR + '/%s.bedintersect.gz')
fulltag_annotation = relfmt(ANNOTATIONS_DIR + '/%s.anno.gz')
fulltag_masked_readids = relfmt(ANNOTATIONS_DIR + '/%s.masked_ids')
fulltag_masked_reads = relfmt(SEQUENCES_DIR + '/full-%s-masked.fa')

shorttag_primary_annotation = relfmt(ANNOTATIONS_DIR + '/short-%s.bedintersect.gz')
shorttag_annotation = relfmt(ANNOTATIONS_DIR + '/short-%s.anno.gz')
shorttag_masked_readids = relfmt(ANNOTATIONS_DIR + '/short-%s.masked_ids')
shorttag_masked_reads = relfmt(SEQUENCES_DIR + '/short-%s-masked.fa')


# =================
# Derived Databases

genomespace_read_database = relfmt(DERIVEDBASES_DIR + '/%s.gspace/done')
genomespace_read_database_dir = relfmt(DERIVEDBASES_DIR + '/%s.gspace')
genomespace_refseq_counts = relfmt(DERIVEDBASES_DIR + '/%s.gspace/refseq.pickle')
genomespace_all_expressed_transcripts = rel(DERIVEDBASES_DIR + '/gspace-expressed-transcripts')

reftranscriptome_sequences = rel(REFTRANSCRIPTOME_DIR + '/reftranscriptome.fa')
reftranscriptome_cds_anno = rel(REFTRANSCRIPTOME_DIR + '/reftranscriptome.bed')
reftranscriptome_gmap_index = rel(REFTRANSCRIPTOME_DIR + '/reftranscriptome/reftranscriptome.genomecomp')
reftranscriptome_dbname = 'reftranscriptome'


# ===========================
# Anti-transcriptome Analyses

fulltag_transcriptome_alignment_sam = relfmt(TRANSCRIPTOMEALN_DIR + '/full-%s.sam.gz')
fulltag_transcriptomic_besthits_sam = relfmt(TRANSCRIPTOMEALN_DIR + '/full-%s-besthits.sam.gz')
fulltag_transcriptome_alignment_gmap = relfmt(TRANSCRIPTOMEALN_DIR + '/full-%s.gmap.gz')
fulltag_transcriptomic_besthits_gmap = relfmt(TRANSCRIPTOMEALN_DIR + '/full-%s-besthits.gmap.gz')
tspace_read_database = relfmt(TRANSCRIPTOMEALN_DIR + '/%s.tspace')


# ==============
# Error Analyses

error_profile_read_level = relfmt(ERRORANALYSIS_DIR + '/%s-readerror.pickle')


# ==========
# Statistics

clip_enrichment_summary = rel(STATISTICS_DIR + '/clip-refseq-enrichment.csv')
total_read_class_stats = relfmt(STATISTICS_DIR + '/classprop.%s.csv')

