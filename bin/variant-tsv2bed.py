#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert the .vcf TSV file to a bed file
"""
import csv
import sys
import argparse

from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)
"""
https://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python
"""

def main(**kwargs):
    """
    Main control function for the script
    """
    input_file = kwargs.pop('input_file', None)
    output_file = kwargs.pop('output_file', None)

    # load input/output file handles
    if input_file:
        fin = open(input_file)
    else:
        fin = sys.stdin

    if output_file:
        fout = open(output_file, "w")
    else:
        fout = sys.stdout

    # start processing input
    reader = csv.DictReader(fin, delimiter = '\t')
    fieldnames = reader.fieldnames
    writer = csv.writer(fout, delimiter = '\t')
    for row in reader:
        chrom = row['CHROM']
        pos = int(row['POS'])
        ref = row['REF']
        alt = row['ALT']

        alt_len = len(alt)
        end = pos + alt_len
        start = pos

        row = [chrom, start, end]
        writer.writerow(row)

    fout.close()
    fin.close()

def parse():
    """
    Parses script args
    """
    parser = argparse.ArgumentParser(description='Filters the ANNOVAR annotation .tsv table for usage with IGV snapshots')
    parser.add_argument("-i", default = None, dest = 'input_file', help="Input file")
    parser.add_argument("-o", default = None, dest = 'output_file', help="Output file")
    args = parser.parse_args()

    main(**vars(args))

if __name__ == '__main__':
    parse()
