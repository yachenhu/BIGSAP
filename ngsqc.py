#!/usr/bin/python3

"""FASTQ quality control."""

import argparse
import sys
import os
import gzip
import subprocess
import random
import math

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '-in_pe', required=True, nargs=2, metavar=('R1.in.fq', 'R2.in.fq'),
        help="Input files of paired-end reads"
    )
    parser.add_argument(
        '-out_pe', required=True, nargs=2, metavar=('R1.out.fq', 'R2.out.fq'),
        help="Output files of paired-end reads"
    )
    parser.add_argument(
        '-max_length', type=int, metavar='INT',
        help=("Reads longer than max_length will be trimmed "
              "to it from 5' end, keep all if not set")
    )
    parser.add_argument(
        '-min_length', type=int, metavar='INT',
        help=("Reads shorter than min_length will be discarded "
              "keep all if not set")
    )
    parser.add_argument(
        '-phred', default=20, type=int, metavar='INT',
        help="Minimum phred score to be HQ [20]"
    )
    parser.add_argument(
        '-percent', default=70, type=float, metavar='FLOAT',
        help="Minimum percent of bases to have [-phred] quality [70]"
    )
    parser.add_argument(
        '-sample', type=int, metavar='INT',
        help="Number of paired reads to sample, keep all if not set"
    )
    parser.add_argument(
        '-compress', action='store_true',
        help="Compress output with GZIP"
    )
    return parser.parse_args()

def gzip_test(fq_file):
    cmd = ['gzip', '-t', fq_file]
    subprocess.run(cmd, stderr=subprocess.PIPE, check=True)

def read(fq_file):
    try:
        gzip_test(fq_file)
    except subprocess.CalledProcessError:
        with open(fq_file, 'r') as fh:
            return list(SeqIO.parse(fh, 'fastq'))
    else:
        with gzip.open(fq_file, 'rt') as fh:
            return list(SeqIO.parse(fh, 'fastq'))

def trim(read, max_length):
    if len(read) > max_length:
        return read[:max_length]
    else:
        return read

def is_high_quality(read, phred, percent):
    scores = read.letter_annotations['phred_quality']
    hq_scores = [n for n in scores if n >= phred]
    if len(hq_scores) / len(scores) * 100 >= percent:
        return True
    else:
        return False

def write(reads, out_file, compress):
    if compress:
        with gzip.open(out_file, 'wt') as fh:
            #SeqIO.write(reads, out_file, 'fastq') not work!
            for read in reads:
                fh.write(read.format('fastq'))
    else:
        with open(out_file, 'w') as fh:
            SeqIO.write(reads, out_file, 'fastq')

def main():
    args = parse_args()
    # Read.
    r1_reads = read(args.in_pe[0])
    r2_reads = read(args.in_pe[1])
    if len(r1_reads) != len(r2_reads):
        return "Error: Number of reads in R1 and R2 is unequal."
    # Begin QC.
    hq_pe_reads = [] # list of read tuple
    for r1_read, r2_read in zip(r1_reads, r2_reads):
        # Trim.
        if args.max_length:
            r1_read = trim(r1_read, args.max_length)
            r2_read = trim(r2_read, args.max_length)
        # Is too short?
        if args.min_length:
            if (len(r1_read) < args.min_length or
                len(r2_read) < args.min_length):
                continue
        # Is low quality?
        if not (is_high_quality(r1_read, args.phred, args.percent) and
                is_high_quality(r2_read, args.phred, args.percent)):
            continue
        hq_pe_reads.append((r1_read, r2_read))
    # Sample.
    if args.sample:
        random.seed(100)
        try:
            hq_pe_reads = random.sample(hq_pe_reads, k=args.sample)
        except ValueError:
            return("Error: Number of paired HQ reads in " +
                   "{} and {} ".format(args.in_pe[0], args.in_pe[1]) +
                   "is less than set value {}".format(args.sample))
    # Write
    for i in range(2):
        out_reads = [r[i] for r in hq_pe_reads]
        write(out_reads, args.out_pe[i], args.compress)

if __name__ == '__main__':
    sys.exit(main())
