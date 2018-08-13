#!/usr/bin/env python

from biotools import bbtools
import multiprocessing
from Bio import SeqIO
from subprocess import PIPE
import subprocess
import argparse
import logging
import shutil  # Not used right now, but will be used for cleanup at some point
import pysam
import os

# TODO: Assemblies with spaces in headers cause problems. Make sure that gets fixed.
# Also, all of this may or may not be doable (and much faster) with regular SNP callers. Look into bcftools/etc
# and see what you can find.


def write_to_logfile(logfile, out, err, cmd):
    with open(logfile, 'a+') as outfile:
        outfile.write('Command used: {}\n\n'.format(cmd))
        outfile.write('STDOUT: {}\n\n'.format(out))
        outfile.write('STDERR: {}\n\n'.format(err))


def run_cmd(cmd):
    p = subprocess.Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    out = out.decode('utf-8')
    err = err.decode('utf-8')
    return out, err


def trim(forward_reads, reverse_reads, outdir):
    """
    Given a set of forward and reverse reads, trims with bbduk and places files called trimmed_R1.fastq.gz and
    trimmed_R2.fastq.gz into outdir - May rename these at some point
    :param forward_reads: Full path to forward reads. Can be gzipped, bzipped, or uncomporessed
    :param reverse_reads: Full path to reverse reads. Can be gzipped, bzipped, or uncomporessed
    :param outdir: Full path to output directory. Must already be craeted.
    """
    forward_out = os.path.join(outdir, 'trimmed_R1.fastq.gz')
    reverse_out = os.path.join(outdir, 'trimmed_R2.fastq.gz')
    cmd = 'bbduk.sh in={forward_in} in2={reverse_in} out={forward_out} out2={reverse_out} qtrim=w trimq=10' \
          ' ref=adapters minlength=50'.format(forward_in=forward_reads,
                                              reverse_in=reverse_reads,
                                              forward_out=forward_out,
                                              reverse_out=reverse_out)
    out, err = run_cmd(cmd)
    write_to_logfile(logfile=os.path.join(outdir, 'log.txt'),
                     out=out,
                     err=err,
                     cmd=cmd)


def create_bam(forward_reads, reverse_reads, reference_fasta, outdir):
    # Map with bbmap default settings. This may or may not work.
    # Getting some real weird results with this - try using bowtie2 and see what happens
    # Bowtie2 gives really similar results with default settings. Will need to play with parameters lots tomorrow to
    # try to figure out what's going on.

    cmd = 'samtools faidx {reference_fasta}'.format(reference_fasta=reference_fasta)
    out, err = run_cmd(cmd)
    write_to_logfile(logfile=os.path.join(outdir, 'log.txt'),
                     out=out,
                     err=err,
                     cmd=cmd)

    bbtools.bbmap(reference=reference_fasta,
                  forward_in=forward_reads,
                  reverse_in=reverse_reads,
                  out_bam=os.path.join(outdir, 'aligned.bam'),
                  subfilter=1)  # Limiting to one substitution allowed per read mapped might help - TBD
    # Also sort and index the bamfile so pysam will be happy with us.
    cmd = 'samtools sort {bamfile} -o {sorted_bamfile}'.format(bamfile=os.path.join(outdir, 'aligned.bam'),
                                                               sorted_bamfile=os.path.join(outdir, 'aligned_sorted.bam'))
    out, err = run_cmd(cmd)
    write_to_logfile(logfile=os.path.join(outdir, 'log.txt'),
                     out=out,
                     err=err,
                     cmd=cmd)

    cmd = 'samtools index {sorted_bamfile}'.format(sorted_bamfile=os.path.join(outdir, 'aligned_sorted.bam'))
    out, err = run_cmd(cmd)
    write_to_logfile(logfile=os.path.join(outdir, 'log.txt'),
                     out=out,
                     err=err,
                     cmd=cmd)


def has_two_high_quality_bases(list_of_scores):
    """
    Finds if a site has at least two bases of very high quality (>35), enough that it can be considered
    fairly safe to say that base is actually there.
    :param list_of_scores: List of quality scores as integers.
    :return: True if site has at least two bases >= 35 phred score, false otherwise
    """
    # TODO: This is a very inelegant way of doing this. Refine.
    quality_bases_count = 0
    for score in list_of_scores:
        if score >= 35:
            quality_bases_count += 1
        if quality_bases_count >= 2:
            return True
    return False


def find_if_multibase(column):
    """
    Finds if a position in a pileup has more than one base present.
    :param column: A pileupColumn generated by pysam
    :return: If position has more than one base, a dictionary with counts for the bases. Otherwise, returns
    empty dictionary
    """
    base_dict = dict()
    # Sometimes the qualities come out to ridiculously high (>70) values. Looks to be because sometimes reads
    # are overlapping and the qualities get summed for overlapping bases. Issue opened on pysam.
    base_qualities = dict()
    for read in column.pileups:
        if read.query_position is not None:  # Not entirely sure why this is sometimes None, but it causes bad stuff
            base = read.alignment.query_sequence[read.query_position]
            quality = read.alignment.query_qualities[read.query_position]
            if base not in base_qualities:
                base_qualities[base] = [quality]
            else:
                base_qualities[base].append(quality)
            if base in base_dict:
                base_dict[base] += 1
            else:
                base_dict[base] = 1
    # Assume that things that are truly indicative of fast evolution at least 20 percent of population?
    # Assume things that only have count of 1 are sequencing errors.
    appropriate_ratio = True
    if len(base_dict) == 1:
        appropriate_ratio = False
    for base in base_dict:
        if base_dict[base]/column.n >= 0.8 or base_dict[base]/column.n <= 0.2 or base_dict[base] == 1:
            appropriate_ratio = False
    if appropriate_ratio is False:
        base_dict = dict()
    # At this point we know that bases are in an appropriate(ish) ratio.
    # Now check that at least two bases for each of the bases present are very high (>35) quality.
    if base_dict:
        high_quality_bases = True
        for base in base_qualities:
            if has_two_high_quality_bases(base_qualities[base]) is False:
                high_quality_bases = False
        if high_quality_bases is False:
            base_dict = dict()
    return base_dict


def get_contig_names(fasta_file):
    """
    Gets contig names from a fasta file using SeqIO.
    :param fasta_file: Full path to uncompressed, fasta-formatted file
    :return: List of contig names.
    """
    contig_names = list()
    for contig in SeqIO.parse(fasta_file, 'fasta'):
        contig_names.append(contig.id)
    return contig_names


def base_dict_to_string(base_dict):
    """
    Converts a dictionary to a string. {'C': 12, 'A':4} gets converted to C:12;A:4
    :param base_dict: Dictionary of bases and counts created by find_if_multibase
    :return: String representing that dictionary.
    """
    outstr = ''
    for base in base_dict:
        outstr += base + ':' + str(base_dict[base]) + ';'
    return outstr[:-1]


def get_average_coverage(bamfile, contig_name):
    # Count coverage gives 4 lists, showing coverage at each position in a contig for A, T, C, and G.
    # To get average coverage, get a list that has sums of the 4 different bases at each position, and then
    # sum that list and find average.
    arrays = bamfile.count_coverage(contig=contig_name)
    list_of_arrays = list()
    for array in arrays:
        list_of_arrays.append(array)
    # Courtesy of: https://stackoverflow.com/questions/14050824/add-sum-of-values-of-two-lists-into-new-list
    summed_list = [sum(x) for x in zip(*list_of_arrays)]
    mean_coverage = sum(summed_list)/len(summed_list)
    return mean_coverage


def read_contig(contig_name, bamfile_name, reference_fasta):
    bamfile = pysam.AlignmentFile(bamfile_name, 'rb')
    multibase_position_dict = dict()
    lines_to_write = list()
    average_coverage = get_average_coverage(bamfile, contig_name)
    # THese parameters seem to be fairly undocumented with pysam, but I think that they should make the output
    # that I'm getting to match up with what I'm seeing in Tablet.
    for column in bamfile.pileup(contig_name,
                                 stepper='samtools', ignore_orphans=False, fastafile=pysam.FastaFile(reference_fasta),
                                 min_base_quality=0):
        base_dict = find_if_multibase(column)
        if base_dict:
            if column.reference_name in multibase_position_dict:
                multibase_position_dict[column.reference_name].append(column.pos)
            else:
                multibase_position_dict[column.reference_name] = [column.pos]
            lines_to_write.append('{},{},{},{},{}\n'.format(column.reference_name,
                                                            column.pos,
                                                            base_dict_to_string(base_dict),
                                                            average_coverage,
                                                            sum(base_dict.values())))
    bamfile.close()
    return lines_to_write


def read_bamfile(sorted_bamfile, outfile, reference_fasta):
    contig_names = get_contig_names(reference_fasta)
    # Make a threads option eventually!
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    bamfile_list = [sorted_bamfile] * len(contig_names)
    reference_fasta_list = [reference_fasta] * len(contig_names)
    # Name this better!
    contig_multi_info = pool.starmap(read_contig, zip(contig_names, bamfile_list, reference_fasta_list), chunksize=1)
    pool.close()
    pool.join()
    with open(outfile, 'w') as f:
        f.write('Contig,Position,Bases,ReadCoverage,BaseCoverage\n')
        for item in contig_multi_info:
            for line in item:
                f.write(line)


def filter_close_snvs(details_file, outfile, snv_window_size, coverage_filter, min_coverage):
    multi_positions = 0
    with open(details_file) as details:
        lines = details.readlines()

    with open(outfile, 'w') as f:
        f.write('Contig,Position,Bases,ReadCoverage,BaseCoverage,Status\n')
        for i in range(1, len(lines)):
            current_line = lines[i].rstrip().split(',')
            current_contig = current_line[0]
            current_position = current_line[1]
            current_bases = current_line[2]
            current_read_coverage = float(current_line[3])
            current_base_coverage = float(current_line[4])
            # Have to check each line against next and previous lines. Special cases for next and previous at end/
            # beginning of file
            if i < len(lines) - 1:
                next_line = lines[i + 1].rstrip().split(',')
                next_contig = next_line[0]
                next_position = next_line[1]
            else:
                next_contig = None
            if i > 1:
                previous_line = lines[i - 1].rstrip().split(',')
                previous_contig = previous_line[0]
                previous_position = previous_line[1]
            else:
                previous_contig = None
            # The case for the vast majority of the file - line is between two others.
            filter_reason = ''
            valid = True
            if next_contig is not None and previous_contig is not None:
                # Check validity: Current postion can't be on same contig and within 500 bases as next or previous
                if current_contig == next_contig and (int(next_position) - int(current_position)) < snv_window_size:
                    valid = False
                    filter_reason = 'SNVDensity'
                elif current_contig == previous_contig and (int(current_position) - int(previous_position)) < snv_window_size:
                    valid = False
                    filter_reason = 'SNVDensity'
            # When we're on first line of file, only need to check if next position, as previous doesn't exist.
            elif previous_contig is None and next_contig is not None:
                if current_contig == next_contig and (int(next_position) - int(current_position)) < snv_window_size:
                    valid = False
                    filter_reason = 'SNVDensity'
            # Final case: last line of file
            elif previous_contig is not None and next_contig is None:
                if current_contig == previous_contig and (int(current_position) - int(previous_position)) < snv_window_size:
                    valid = False
                    filter_reason = 'SNVDensity'

            # Also: check that base coverage isn't significantly (1.5X by default) more than average coverage.
            # Can indicate that two (or more) slightly different genes have been assembled into one by SPAdes
            if current_base_coverage/current_read_coverage > coverage_filter:
                valid = False
                filter_reason = 'CoverageHigh'

            if current_base_coverage < min_coverage:
                valid = False
                filter_reason = 'CoverageLow'

            if valid:
                multi_positions += 1
                f.write('{},{},{},{},{},{}\n'.format(current_contig,
                                                     current_position,
                                                     current_bases,
                                                     current_read_coverage,
                                                     current_base_coverage,
                                                     'Valid'))
            else:
                f.write('{},{},{},{},{},{}\n'.format(current_contig,
                                                     current_position,
                                                     current_bases,
                                                     current_read_coverage,
                                                     current_base_coverage,
                                                     'Filtered-' + filter_reason))

    return multi_positions


def cleanup(directory):
    os.system('rm {directory}/*.fastq.gz {directory}/*bam*'.format(directory=directory))


def main():
    logging.basicConfig(format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
                        level=logging.INFO,
                        datefmt='%Y-%m-%d %H:%M:%S')
    parser = argparse.ArgumentParser(description='Attempts to find evolutionary rate of a strain, when provided'
                                                 ' with paired-end FASTQ files and an assembly for that strain.')
    parser.add_argument('-a', '--assembly',
                        type=str,
                        required=True,
                        help='Full path to your FASTA-formatted assembly. Only works with SPAdes assemblies right now.'
                             ' SKESA causes breakpoints at most sites we\'re looking for. Other assemblers have'
                             ' not been tested.')
    parser.add_argument('-1', '--forward_reads',
                        type=str,
                        required=True,
                        help='Full path to your forward (_R1) read file.')
    parser.add_argument('-2', '--reverse_reads',
                        type=str,
                        required=True,
                        help='Full path to your reverse(_R2) read file.')
    parser.add_argument('-o', '--outdir',
                        type=str,
                        required=True,
                        help='Output directory. Must not currently exist.')
    parser.add_argument('-m', '--min_coverage',
                        default=10,
                        type=int,
                        help='Minimum amount of coverage a site must have before it will be considered for analysis. '
                             'Defaults to 10.')
    parser.add_argument('-w', '--filter_density_window',
                        default=500,
                        type=int,
                        help='Window size for identifying high-density SNV regions that will be discarded from '
                             'analysis. Defaults to 500, meaning any region of 500 or fewer base pairs with more than '
                             'one SNV will not be considered.')
    parser.add_argument('-c', '--coverage_filter',
                        type=float,
                        default=1.5,
                        help='Filter for coverage - sites that look heterogenous but have higher than average '
                             'coverage usually result from two genes that are very close to identical that get '
                             'merged into one by SPAdes. This filter stops that. Defaults to discarding if coverage '
                             'is more than 1.5x the average for the contig.')
    args = parser.parse_args()

    if os.path.isdir(args.outdir):
        logging.error('ERROR: Output directory specified already exists! Must be a new directory. Try again.')
        quit(code=1)
    else:
        os.makedirs(args.outdir)
    # Steps in this analysis:
    # 1) Trim forward and reverse reads. I think this is a good idea, but maybe we'll end up leaving this out.

    # This gets us trimmed_R1.fastq.gz and trimmed_R2.fastq.gz in our temporary directory.
    logging.info('Trimming reads...')
    trim(forward_reads=args.forward_reads,
         reverse_reads=args.reverse_reads,
         outdir=args.outdir)
    # 2) Align trimmed reads back to reference. Start off using bbmap for this, but maybe try bowtie2 as well.

    logging.info('Aligning reads back to reference...')
    # We get a sorted and indexed bamfile called aligned_sorted.bam in our temporary directory.
    create_bam(forward_reads=os.path.join(args.outdir, 'trimmed_R1.fastq.gz'),
               reverse_reads=os.path.join(args.outdir, 'trimmed_R2.fastq.gz'),
               reference_fasta=args.assembly,
               outdir=args.outdir)
    # 3) Parse through each position of the BAM file and find sites that are heterogenous-ish, indicating the colony
    # was evolving quickly.
    logging.info('Parsing BAM file...')
    read_bamfile(sorted_bamfile=os.path.join(args.outdir, 'aligned_sorted.bam'),
                 outfile=os.path.join(args.outdir, 'details.csv'),
                 reference_fasta=args.assembly)

    # This gets us a list of sites that are appropriately heterogenous (pending me figuring out what cutoffs should
    # actually be getting used). A lot of these sites tend to cluster together - my interpretation is that this is
    # because there are potentially paralogous genes present, so a lot of reads from gene1 map to gene2, making it
    # look like it's evolved really quickly, but isn't actually the case.

    # To resolve this: Use some sort of sliding window to find these high-density regions. Use SNVPhyl settings
    # (only 1SNV/500bp) to start.
    logging.info('Filtering...')
    multi_positions = filter_close_snvs(details_file=os.path.join(args.outdir, 'details.csv'),
                                        outfile=os.path.join(args.outdir, 'details_filtered.csv'),
                                        snv_window_size=args.filter_density_window,
                                        coverage_filter=args.coverage_filter,
                                        min_coverage=args.min_coverage)
    cleanup(args.outdir)
    logging.info('Complete! Total multi positions found: {}'.format(multi_positions))


if __name__ == '__main__':
    main()
