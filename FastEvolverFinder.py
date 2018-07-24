#!/usr/bin/env python

from biotools import bbtools
import subprocess
import argparse
import shutil
import pysam
import os


def trim(forward_reads, reverse_reads, outdir):
    forward_out = os.path.join(outdir, 'trimmed_R1.fastq.gz')
    reverse_out = os.path.join(outdir, 'trimmed_R2.fastq.gz')
    cmd = 'bbduk.sh in={forward_in} in2={reverse_in} out={forward_out} out2={reverse_out} qtrim=w trimq=10' \
          ' ref=adapters minlength=50'.format(forward_in=forward_reads,
                                              reverse_in=reverse_reads,
                                              forward_out=forward_out,
                                              reverse_out=reverse_out)
    subprocess.call(cmd, shell=True)


def create_bam(forward_reads, reverse_reads, reference_fasta, outdir):
    # Map with bbmap default settings. This may or may not work.
    bbtools.bbmap(reference=reference_fasta,
                  forward_in=forward_reads,
                  reverse_in=reverse_reads,
                  out_bam=os.path.join(outdir, 'aligned.bam'))
    # Also sort and index the bamfile so pysam will be happy with us.
    cmd = 'samtools sort {bamfile} -o {sorted_bamfile}'.format(bamfile=os.path.join(outdir, 'aligned.bam'),
                                                               sorted_bamfile=os.path.join(outdir, 'aligned_sorted.bam'))
    subprocess.call(cmd, shell=True)
    cmd = 'samtools index {sorted_bamfile}'.format(sorted_bamfile=os.path.join(outdir, 'aligned_sorted.bam'))
    subprocess.call(cmd, shell=True)


def find_if_multibase(column):
    # Try with 3 percent cutoff, only 2 bases present?
    base_dict = dict()
    for read in column.pileups:
        if read.query_position is not None:  # Not entirely sure why this is sometimes None, but it causes bad stuff
            base = read.alignment.query_sequence[read.query_position]
            if base in base_dict:
                base_dict[base] += 1
            else:
                base_dict[base] = 1
    # Assume that things that are truly indicative of fast evolution at least 10 percent of population?
    # Assume things that only have count of 1 are sequencing errors. TODO: See if we can get quality to inform this.
    appropriate_ratio = True
    if len(base_dict) == 1:
        appropriate_ratio = False
    for base in base_dict:
        if base_dict[base]/column.n >= 0.9 or base_dict[base]/column.n <= 0.1 or base_dict[base] == 1:
            appropriate_ratio = False
    # if appropriate_ratio:
    #     print(base_dict, end=' ')
    #     print(column.pos, end=' ')
    #     print(column.reference_name)
    if appropriate_ratio is False:
        base_dict = dict()
    return base_dict


def base_dict_to_string(base_dict):
    outstr = ''
    for base in base_dict:
        outstr += base + ':' + str(base_dict[base]) + ';'
    return outstr[:-1]


def read_bamfile(sorted_bamfile, outfile):
    multi_base_positions = 0
    bamfile = pysam.Samfile(sorted_bamfile, 'rb')
    # TODO: Figure out if the multiprocessing is actually possible. Currently don't think it is due to some
    # limitation with pysam
    # pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    # multibase = pool.imap(find_if_multibase, bamfile.pileup(), chunksize=1000)
    # pool.close()
    # pool.join()
    with open(outfile, 'w') as f:
        f.write('Contig,Position,Bases\n')
        multibase_position_dict = dict()
        for column in bamfile.pileup():
            base_dict = find_if_multibase(column)
            if base_dict:
                multi_base_positions += 1
                if column.reference_name in multibase_position_dict:
                    multibase_position_dict[column.reference_name].append(column.pos)
                else:
                    multibase_position_dict[column.reference_name] = [column.pos]
                f.write('{},{},{}\n'.format(column.reference_name, column.pos, base_dict_to_string(base_dict)))
            # if total_genome_length % 200000 == 0:  # Give me some sporadic updates.
            #     print('Genome length ' + str(total_genome_length))
            #     print('Multi base positions ' + str(multi_base_positions))
    bamfile.close()
    # print('Genome length ' + str(total_genome_length))
    # print('Multi base positions ' + str(multi_base_positions))
    # print(multibase_position_dict)
    return multibase_position_dict


def filter_close_snvs(details_file, outfile):
    multi_positions = 0
    with open(details_file) as details:
        lines = details.readlines()

    with open(outfile, 'w') as f:
        f.write('Contig,Position,Bases,Status\n')
        for i in range(1, len(lines)):
            current_line = lines[i].rstrip().split(',')
            current_contig = current_line[0]
            current_position = current_line[1]
            current_bases = current_line[2]
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
            valid = True
            if next_contig is not None and previous_contig is not None:
                # Check validity: Current postion can't be on same contig and within 500 bases as next or previous
                if current_contig == next_contig and (int(next_position) - int(current_position)) < 500:
                    valid = False
                elif current_contig == previous_contig and (int(current_position) - int(previous_position)) < 500:
                    valid = False
            # When we're on first line of file, only need to check if next position, as previous doesn't exist.
            elif previous_contig is None and next_contig is not None:
                if current_contig == next_contig and (int(next_position) - int(current_position)) < 500:
                    valid = False
            # Final case: last line of file
            elif previous_contig is not None and next_contig is None:
                if current_contig == previous_contig and (int(current_position) - int(previous_position)) < 500:
                    valid = False

            if valid:
                multi_positions += 1
                f.write('{},{},{},{}\n'.format(current_contig, current_position, current_bases, 'Valid'))
            else:
                f.write('{},{},{},{}\n'.format(current_contig, current_position, current_bases, 'Filtered'))

    return multi_positions


def main():
    parser = argparse.ArgumentParser(description='Attempts to find evolutionary rate of a strain, when provided'
                                                 ' with paired-end FASTQ files and an assembly for that strain.')
    parser.add_argument('-a', '--assembly',
                        type=str,
                        required=True,
                        help='Full path to your FASTA-formatted assembly.')
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
    args = parser.parse_args()

    if os.path.isdir(args.outdir):
        print('ERROR: Output directory specified already exists! Must be a new directory. Try again.')
        # quit(code=1)
    else:
        os.makedirs(args.outdir)
    # Steps in this analysis:
    # 1) Trim forward and reverse reads. I think this is a good idea, but maybe we'll end up leaving this out.

    # This gets us trimmed_R1.fastq.gz and trimmed_R2.fastq.gz in our temporary directory.
    trim(forward_reads=args.forward_reads,
         reverse_reads=args.reverse_reads,
         outdir=args.outdir)
    # 2) Align trimmed reads back to reference. Start off using bbmap for this, but maybe try bowtie2 as well.

    # We get a sorted and indexed bamfile called aligned_sorted.bam in our temporary directory.
    create_bam(forward_reads=os.path.join(args.outdir, 'trimmed_R1.fastq.gz'),
               reverse_reads=os.path.join(args.outdir, 'trimmed_R2.fastq.gz'),
               reference_fasta=args.assembly,
               outdir=args.outdir)
    # 3) Parse through each position of the BAM file and find sites that are heterogenous-ish, indicating the colony
    # was evolving quickly.
    multi_position_dict = read_bamfile(sorted_bamfile=os.path.join(args.outdir, 'aligned_sorted.bam'),
                                       outfile=os.path.join(args.outdir, 'details.csv'))

    # This gets us a list of sites that are appropriately heterogenous (pending me figuring out what cutoffs should
    # actually be getting used). A lot of these sites tend to cluster together - my interpretation is that this is
    # because there are potentially paralogous genes present, so a lot of reads from gene1 map to gene2, making it
    # look like it's evolved really quickly, but isn't actually the case.

    # To resolve this: Use some sort of sliding window to find these high-density regions. Use SNVPhyl settings
    # (only 1SNV/500bp) to start. TODO: Maybe implement some BLAST search to see if regions found are actually
    multi_positions = filter_close_snvs(details_file=os.path.join(args.outdir, 'details.csv'),
                                        outfile=os.path.join(args.outdir, 'details_filtered.csv'))

    print('Total multi positions found: {}'.format(multi_positions))


if __name__ == '__main__':
    main()
