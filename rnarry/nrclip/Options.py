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
NUM_THREADS = 24
MAX_PARALLEL_DOWNLOADING = 4

GENOME = 'mm9'
SPECIES = 'Mus musculus'
SPECIES3 = 'mmu'
CLIP_SAMPLES = ['CLIP-35L33G', 'CLIP-2J3', 'CLIP-46020']
CLIPCTL_SAMPLES = ['PolyA-None-rep1']
RPF_SAMPLES = ['RPF-Luc-rep1', 'RPF-Lin28a-rep1']
RPFCTL_SAMPLES = ['PolyA-Luc-rep1', 'PolyA-Lin28a-rep1']
ALL_SAMPLES = CLIP_SAMPLES + CLIPCTL_SAMPLES + RPF_SAMPLES + RPFCTL_SAMPLES
SHORTTAG_SAMPLES = RPF_SAMPLES + RPFCTL_SAMPLES

modban_illumina = 'CTGTAGGCACCATCAATTCGTATGCCGTCTTCTGCTTG'
illumina_SRA15 = 'ATCTCGTATGCCGTCTTCTGCTTG'
ADAPTER_SEQ = {
    'CLIP-35L33G': modban_illumina,
    'CLIP-2J3': modban_illumina,
    'CLIP-46020': modban_illumina,
    'PolyA-None-rep1': modban_illumina,
    'RPF-Luc-rep1': illumina_SRA15,
    'RPF-Lin28a-rep1': illumina_SRA15,
    'PolyA-Luc-rep1': illumina_SRA15,
    'PolyA-Lin28a-rep1': illumina_SRA15,
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

GENOMEALN_POSTPROC_ALLOWED_MISMATCHES = 2

# Annotation settings
UCSC_REPEATMASKER_CHROMOSOMES = """
1 10 11 12 13 13_random 14 15 16 16_random 17 17_random 18
7 19 1_random 2 3 3_random 4 4_random 5 5_random 6 7 7_random 8 8_random 9
8 9_random M Un_random X X_random Y Y_random""".split()
# These RNAs are annotated by gtRNAdb and/or rfam.
REPEATMASKER_IGNORE_CLASSES = 'tRNA snRNA scRNA srpRNA'.split()

