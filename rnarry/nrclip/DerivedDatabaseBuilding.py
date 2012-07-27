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


@files([Paths.nr_refseq_list] +
       [Paths.genomespace_refseq_counts(sample)
        for sample in Paths.ALL_SAMPLES],
       Paths.genomespace_all_expressed_transcripts)
@follows(quantitate_refseq_in_gspace)
def make_list_of_expressed_transcripts(inputfiles, outputfile):
    inputlist = ' '.join(inputfiles)

    runproc("""
        $ENV MINDEPTH=$GSPACE_STATS_MINIMUM_RAW_READS \
        $REFSEQCNT_PICK_EXPRESSED $inputlist > $outputfile""", outputfile)


@files([Paths.nr_refseq_db, Paths.genomespace_all_expressed_transcripts,
        Paths.genome_fasta,
        Paths.genomespace_read_database(Options.SNPREF_SAMPLE)],
       [Paths.reftranscriptome_sequences, Paths.reftranscriptome_cds_anno])
@follows(make_list_of_expressed_transcripts)
def generate_SNPfixed_transcriptome(inputfiles, outputfiles):
    nrdb, expressed, genomeseq, gspace = inputfiles
    outseq, outanno = outputfiles

    runproc("""
        $GENFASTA_MUTATED_TRANSCRIPTOME $nrdb $expressed $genomeseq $gspace \
            $outseq $outanno""", outputfiles)
    runproc("$SAMTOOLS faidx $outseq", outputfiles)
    runproc("$FASIZE -detailed $outseq > $outseq.size", outputfiles)


@files(Paths.reftranscriptome_sequences, Paths.reftranscriptome_gmap_index)
@follows(generate_SNPfixed_transcriptome)
def generate_gsnap_transcriptome_index(inputfile, outcheck):
    runproc('$GMAP_BUILD -d $reftranscriptome_dbname -D $REFTRANSCRIPTOME_DIR '
            '-k 12 $inputfile', outcheck)


def tasks():
    return [
        build_genomespace_read_database,
        quantitate_refseq_in_gspace,
        make_list_of_expressed_transcripts,
        generate_SNPfixed_transcriptome,
        generate_gsnap_transcriptome_index,
    ]
