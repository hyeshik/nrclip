#
# rnarry.nrclip.ErrorAnalysis
#  - Error profile analysis for CLIP-seq samples
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
from rnarry.nrclip import Paths, Options, TranscriptomeAnalysis
from rnarry.nrclip.PipelineControl import *
from rnarry.sequtils import get_first_sequence_length
from itertools import chain
import os

FASTQ_FORMAT = {
    'illumina1.5': 'fastq-illumina',
    'illumina1.8': 'fastq-sanger',
}

@files(for_each([Paths.fulltag_transcriptomic_besthits_gmap,
                 Paths.original_sequence_reads],
                Paths.error_profile_read_level,
                Paths.ALLCLIP_SAMPLES))
@follows(TranscriptomeAnalysis.resolve_transcriptomic_multihits_gmap)
def count_alignment_errors(inputfiles, outputfile, sample):
    mapping, origseqs = inputfiles
    format = FASTQ_FORMAT[Options.QUALITY_SCALE[sample]]
    readlength = get_first_sequence_length(origseqs, format=format)

    runproc('$GMAP_ERROR_PROFILE $mapping $outputfile $readlength',
            outputfile)


@files(for_each(Paths.error_profile_read_level,
                Paths.error_profile_summarized,
                Paths.ALLCLIP_SAMPLES))
@follows(count_alignment_errors)
def summarize_error_profile(inputfile, outputfile, sample):
    runproc('$SUMMARIZE_ERROR_PROFILE $inputfile $outputfile', outputfile)


@files(for_each([Paths.reftranscriptome_sequences,
                 Paths.fulltag_transcriptomic_besthits_gmap,
                 Paths.error_profile_read_level],
                Paths.clipsim_input_data_pack,
                Paths.ALLCLIP_SAMPLES))
@follows(TranscriptomeAnalysis.resolve_transcriptomic_multihits_gmap)
@follows(count_alignment_errors)
def prepare_clip_sim_input(inputfiles, outputfile, sample):
    inputs = ' '.join(inputfiles)
    runproc('$CLIPSIM_PREPARE_INPUTS $inputs $outputfile', outputfile)


@files(for_each(Paths.clipsim_input_data_pack,
                [Paths.clipsim_permutated_del,
                 Paths.clipsim_permutated_mod,
                 Paths.clipsim_permutated_moddel,
                 Paths.clipsim_permutated_entropy],
                Paths.ALLCLIP_SAMPLES))
@follows(prepare_clip_sim_input)
@jobs_limit(1, 'exclusive')
def permutate_clip_alignments(inputfile, outputfiles, sample):
    outputprefix = os.path.commonprefix(outputfiles)
    runproc('$CROSSFEST -i $inputfile -o $outputprefix -t $NUM_THREADS '
            '-d $CROSSFEST_MINIMUM_COVERAGE', outputfiles)


@files(for_each([Paths.fulltag_transcriptomic_besthits_gmap,
                 Paths.reftranscriptome_sequences],
                [Paths.clipsim_real_nonzero_positions,
                 Paths.clipsim_real_del_scores,
                 Paths.clipsim_real_mod_scores,
                 Paths.clipsim_real_moddel_scores,
                 Paths.clipsim_real_entropy_scores],
                Paths.ALLCLIP_SAMPLES))
@follows(TranscriptomeAnalysis.resolve_transcriptomic_multihits_gmap)
def calculate_clip_scores_real(inputfiles, outputfiles, sample):
    inputs = ' '.join(inputfiles)
    outputprefix = os.path.commonprefix(outputfiles)
    runproc('$CLIPSIM_REALDATA_DISTS $inputs $outputprefix', outputfiles)


@files(list(chain(*[for_each([Paths.clipsim_files_by_method[method][0],
                              Paths.clipsim_files_by_method[method][1]],
                             Paths.clipsim_files_by_method[method][2],
                             Paths.ALLCLIP_SAMPLES)
                    for method in Paths.MUTATION_RATE_METHODS])))
@follows(permutate_clip_alignments)
@follows(calculate_clip_scores_real)
def calculate_cutoffs_by_fdr(inputfiles, outputfile, sample):
    inputprefix = os.path.commonprefix(inputfiles)
    runproc('$CLIPSIM_CALC_FDR_CURVE $inputprefix $outputfile', outputfile)


def tasks():
    return [
        count_alignment_errors,
        summarize_error_profile,
        prepare_clip_sim_input,
        permutate_clip_alignments,
        calculate_clip_scores_real,
        calculate_cutoffs_by_fdr,
    ]

