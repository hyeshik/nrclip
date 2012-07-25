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
from rnarry.nrclip import Paths, Options
from rnarry.nrclip.PipelineControl import *


@files(Paths.early_filter_fasta, Paths.early_filter_index_check)
def generate_gsnap_early_filter_index(inputfile, outcheck):
    runproc('$GMAP_BUILD -d $early_filter_prefix '
            '-D $EXTERNAL_DIR $inputfile', outcheck)


@files(None, Paths.genome_fasta_zipped)
@jobs_limit(Options.MAX_PARALLEL_DOWNLOADING, 'download')
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
@jobs_limit(Options.MAX_PARALLEL_DOWNLOADING, 'download')
def download_refgene(inputfile, outputfile):
    runproc('$WGET -O $outputfile $REFGENE_URL', outputfile)


@files(None, Paths.knowngene_ucsc)
@jobs_limit(Options.MAX_PARALLEL_DOWNLOADING, 'download')
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


@files(None, Paths.mirbase_catalog)
@jobs_limit(Options.MAX_PARALLEL_DOWNLOADING, 'download')
def download_mirbase_catalog(inputfile, outputfile):
    import urllib, re, gzip

    accession_pat = re.compile('ACC.*ID="([^"]*)";')

    with DeleteOnError(outputfile, gzip.open) as output:
        for line in urllib.urlopen(Paths.MIRBASE_URL):
            if line.startswith('#'):
                continue
            fields = line[:-1].split('\t')
            name = accession_pat.findall(fields[8])[0]
            print >> output, '\t'.join([
                'chr'+fields[0], str(int(fields[3])-1), fields[4],
                'miRNA|%s|%s' % (name, name), '.', fields[6]
            ])


@files([(None, Paths.repeatmasker_original(chrom), chrom)
        for chrom in Options.UCSC_REPEATMASKER_CHROMOSOMES])
@jobs_limit(Options.MAX_PARALLEL_DOWNLOADING, 'download')
def download_repeatmasker_catalogs(inputfile, outputfile, chrom):
    import urllib

    url = Paths.REPEATMASKER_URL % (Options.GENOME, chrom)
    print 'Downloading %s ...' % url
    urllib.urlretrieve(url, outputfile)


@files([Paths.repeatmasker_original(chrom)
        for chrom in Options.UCSC_REPEATMASKER_CHROMOSOMES],
       Paths.repeatmasker_catalog)
@follows(download_repeatmasker_catalogs)
def merge_repeatmasker_catalogs(inputfiles, outputfile):
    import gzip, itertools

    igncls = set(Options.REPEATMASKER_IGNORE_CLASSES)
    with DeleteOnError(outputfile, gzip.open) as output:
        for line in itertools.chain(*[gzip.open(inp) for inp in inputfiles]):
            fields = line[:-1].split('\t')
            if fields[12] in igncls:
                continue

            print >> output, '%s\t%s\t%s\t%s|%s|%s\t%s\t%s' % (
                fields[5], fields[6], fields[7],
                fields[11], fields[10], fields[12] + '/' + fields[10],
                fields[1], fields[9])


@files(Paths.refgene_ucsc, Paths.refseq_catalog)
@follows(download_refgene)
def prepare_refseq_catalog(inputfile, outputfile):
    outputdir = os.path.dirname(outputfile)
    runproc('$BUILD_REFSEQ_INDEX $outputdir $inputfile', outputfile)


@files(Paths.genome_fasta, Paths.genome_twobit)
@follows(extract_genome_sequence)
def prepare_genome_2bit(inputfile, outputfile):
    runproc('$FA2TWOBIT $inputfile $outputfile', outputfile)


@files(None, [Paths.rfam_fasta, Paths.rfam_original])
@jobs_limit(Options.MAX_PARALLEL_DOWNLOADING, 'download')
def download_rfam_files(inputfile, outputfiles):
    import urllib

    origurls = dict((os.path.basename(url), url)
                   for url in (Paths.RFAM_FULL_URL, Paths.RFAM_SEQUENCE_URL))

    for outputfile in outputfiles:
        url = origurls[os.path.basename(outputfile)]
        print "Downloading %s ..." % url
        urllib.urlretrieve(url, outputfile)


@files([Paths.rfam_fasta, Paths.rfam_original], Paths.rfam_catalog)
@follows(download_rfam_files)
@follows(prepare_genome_2bit)
def prepare_rfam_catalog(inputfiles, outputfile):
    rfamseq, rfamfull = inputfiles

    runproc("""$BUILD_RFAM_INDEX $TMP_DIR $rfam_fasta $rfam_original \
               "$SPECIES" $genome_twobit $outputfile""", outputfile)


@files(None, Paths.trna_catalog)
@jobs_limit(Options.MAX_PARALLEL_DOWNLOADING, 'download')
def download_trna_catalog(inputfile, outputfile):
    import urllib, re, gzip, StringIO

    inpbuf = StringIO.StringIO(urllib.urlopen(Paths.GTRNADB_URL).read())

    with DeleteOnError(outputfile, gzip.open) as output:
        for line in gzip.GzipFile(fileobj=inpbuf):
            fields = line[:-1].split('\t')
            print >> output, '\t'.join([
                fields[1], fields[2], fields[3],
                'tRNA|tRNA-%s%s|tRNA-%s%s' % (
                    fields[7], fields[8], fields[7], fields[8]),
                fields[5], fields[6]])


@files([
    Paths.mirbase_catalog, Paths.refseq_catalog, Paths.repeatmasker_catalog,
    Paths.rfam_catalog, Paths.trna_catalog], Paths.compiled_catalog)
@follows(download_mirbase_catalog)
@follows(prepare_refseq_catalog)
@follows(merge_repeatmasker_catalogs)
@follows(prepare_rfam_catalog)
@follows(download_trna_catalog)
def compile_catalogs(inputfiles, outputfile):
    inputs = ' '.join(inputfiles)
    runproc("""
        $ZCAT $inputs | sort -k1,1 -k2,3n -k4,4 |
        $GZIP -c - > $outputfile""", outputfile)


def tasks():
    return [
        generate_gsnap_early_filter_index,
        download_zipped_genome_sequence,
        extract_genome_sequence,
        generate_gsnap_genome_index,
        download_refgene,
        download_knowngene,
        build_splicesites_index,
        download_mirbase_catalog,
        download_repeatmasker_catalogs,
        merge_repeatmasker_catalogs,
        prepare_refseq_catalog,
        prepare_genome_2bit,
        download_rfam_files,
        prepare_rfam_catalog,
        download_trna_catalog,
        compile_catalogs,
    ]
