#
# rnarry.nrclip.Options
#  - Options and parameters
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


# ================================
# Project-specific configurations

NUM_PARALLEL = 8
NUM_THREADS = 18
#NUM_THREADS = 8
MAX_PARALLEL_DOWNLOADING = 4

GENOME = 'mm9'
SPECIES = 'Mus musculus'
SPECIES3 = 'mmu'
CLIP_SAMPLES = ['CLIP-35L33G', 'CLIP-2J3', 'CLIP-46020']
CLIPCTL_SAMPLES = ['PolyA-None-rep1', 'PolyA-None-rep2', 'PolyA-None-rep3']
RPF_SAMPLES = ['RPF-Luc-rep1', 'RPF-Luc-rep2', 'RPF-GFP-rep1',
               'RPF-Lin28a-rep1', 'RPF-Lin28a-rep2', 'RPF-Lin28a-rep3']
RPFCTL_SAMPLES = ['PolyA-Luc-rep1', 'PolyA-Luc-rep2', 'PolyA-GFP-rep1',
                  'PolyA-Lin28a-rep1', 'PolyA-Lin28a-rep2', 'PolyA-Lin28a-rep3']
ALL_SAMPLES = CLIP_SAMPLES + CLIPCTL_SAMPLES + RPF_SAMPLES + RPFCTL_SAMPLES
SHORTTAG_SAMPLES = RPF_SAMPLES + RPFCTL_SAMPLES
FULLTAG_SAMPLES = CLIP_SAMPLES + CLIPCTL_SAMPLES
SNPREF_SAMPLE = CLIPCTL_SAMPLES[0]
ALLCLIP_SAMPLES = CLIP_SAMPLES + CLIPCTL_SAMPLES
RPF_PAIRS = [
    ('Luc-r1', ('RPF-Luc-rep1', 'PolyA-Luc-rep1')),
    ('Luc-r2', ('RPF-Luc-rep2', 'PolyA-Luc-rep2')),
    ('GFP-r1', ('RPF-GFP-rep1', 'PolyA-GFP-rep1')),
    ('Lin28a-r1', ('RPF-Lin28a-rep1', 'PolyA-Lin28a-rep1')),
    ('Lin28a-r2', ('RPF-Lin28a-rep2', 'PolyA-Lin28a-rep2')),
    ('Lin28a-r3', ('RPF-Lin28a-rep3', 'PolyA-Lin28a-rep3')),
]

modban_illumina = 'CTGTAGGCACCATCAATTCGTATGCCGTCTTCTGCTTG'
illumina_SRA15 = 'ATCTCGTATGCCGTCTTCTGCTTG'
illumina_TruSeq = 'TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC'
ADAPTER_SEQ = {
    'CLIP-35L33G': modban_illumina,
    'CLIP-2J3': modban_illumina,
    'CLIP-46020': modban_illumina,
    'PolyA-None-rep1': modban_illumina,
    'PolyA-None-rep2': illumina_TruSeq,
    'PolyA-None-rep3': illumina_TruSeq,
    'RPF-Luc-rep1': illumina_SRA15,
    'RPF-Luc-rep2': illumina_TruSeq,
    'RPF-GFP-rep1': illumina_TruSeq,
    'RPF-Lin28a-rep1': illumina_SRA15, # siLin28a-8
    'RPF-Lin28a-rep2': illumina_TruSeq, # siLin28a-8
    'RPF-Lin28a-rep3': illumina_TruSeq, # siLin28a-9
    'PolyA-Luc-rep1': illumina_SRA15,
    'PolyA-Luc-rep2': illumina_TruSeq,
    'PolyA-GFP-rep1': illumina_TruSeq,
    'PolyA-Lin28a-rep1': illumina_SRA15, # siLin28a-8
    'PolyA-Lin28a-rep2': illumina_TruSeq, # siLin28a-8
    'PolyA-Lin28a-rep3': illumina_TruSeq, # siLin28a-9
}
QUALITY_SCALE = {
    'CLIP-35L33G': 'illumina1.5',
    'CLIP-2J3': 'illumina1.5',
    'CLIP-46020': 'illumina1.5',
    'PolyA-None-rep1': 'illumina1.5',
    'PolyA-None-rep2': 'illumina1.8',
    'PolyA-None-rep3': 'illumina1.8',
    'RPF-Luc-rep1': 'illumina1.5',
    'RPF-Luc-rep2': 'illumina1.8',
    'RPF-GFP-rep1': 'illumina1.8',
    'RPF-Lin28a-rep1': 'illumina1.5', # siLin28a-8
    'RPF-Lin28a-rep2': 'illumina1.8', # siLin28a-8
    'RPF-Lin28a-rep3': 'illumina1.8', # siLin28a-9
    'PolyA-Luc-rep1': 'illumina1.5',
    'PolyA-Luc-rep2': 'illumina1.8',
    'PolyA-GFP-rep1': 'illumina1.8',
    'PolyA-Lin28a-rep1': 'illumina1.5', # siLin28a-8
    'PolyA-Lin28a-rep2': 'illumina1.8', # siLin28a-8
    'PolyA-Lin28a-rep3': 'illumina1.8', # siLin28a-9
}

# fastx toolkit settings
FULLTAG_MIN_LENGTH = 20
FULLTAG_MIN_QUALITY = 25
FULLTAG_MIN_QUALITY_PERCENT = 90
SHORTTAG_LENGTH = 27
SHORTTAG_MIN_QUALITY = 25
SHORTTAG_MIN_QUALITY_PERCENT = 90

# GSNAP alignment options
FULLTAG_PREALN_MISMATCHES = '0.1'
SHORTTAG_PREALN_MISMATCHES = '0.05'
FULLTAG_GENOME_MISMATCHES = '0.1'
SHORTTAG_GENOME_MISMATCHES = '0.05'
FULLTAG_TRANSCRIPTOME_MISMATCHES = '0.05' # SNPs are fixed already
SHORTTAG_TRANSCRIPTOME_MISMATCHES = '0.05'

FULLTAG_POSTPROC_ALLOWED_MISMATCHES = 2 # for CLIP-seq
SHORTTAG_POSTPROC_ALLOWED_MISMATCHES = 1 # for RPF

FULLTAG_GMAP_MAXIMUM_INDEL = 3 # for transcriptome analyses

# Annotation settings
UCSC_REPEATMASKER_CHROMOSOMES = """
1 10 11 12 13 13_random 14 15 16 16_random 17 17_random 18
7 19 1_random 2 3 3_random 4 4_random 5 5_random 6 7 7_random 8 8_random 9
8 9_random M Un_random X X_random Y Y_random""".split()
# These RNAs are annotated by gtRNAdb and/or rfam.
REPEATMASKER_IGNORE_CLASSES = 'tRNA snRNA scRNA srpRNA'.split()

# Expression quantification settings
GSPACE_STATS_MINIMUM_RAW_READS = 10

# CLIP options
CROSSFEST_PERMUTATIONS = 1000

