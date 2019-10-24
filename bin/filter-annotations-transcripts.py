#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Filter annotation table to keep only desired transcripts from provided transcript list
"""
import csv
import sys
import argparse

transcript_colname = "Transcript"
NA_char = "."
invalid_chars = ["UNKNOWN", NA_char]
def main(**kwargs):
    """
    Main control function for the script
    """
    input_file = kwargs.pop('input_file', None)
    output_file = kwargs.pop('output_file', None)
    transcripts_file = kwargs.pop('transcripts_file', "transcripts.txt")

    #  read all transcripts
    with open(transcripts_file) as f:
        transcripts = set([ t.strip() for t in f.readlines() ])

    transcripts.add(NA_char)

    # open input/output filehandles
    if input_file:
        fin = open(input_file)
    else:
        fin = sys.stdin

    if output_file:
        fout = open(output_file, "w")
    else:
        fout = sys.stdout

    # start tsv parsing
    reader = csv.DictReader(fin, delimiter = '\t')
    fieldnames = [f for f in reader.fieldnames]
    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = fieldnames)
    writer.writeheader()

    for row in reader:
        if row[transcript_colname] in transcripts:
            writer.writerow(row)

    fout.close()
    fin.close()


def parse():
    """
    Parses script args
    """
    parser = argparse.ArgumentParser(description='Filters the annotation table to output only annotations that match desired transcripts. Keeps entries with NA or empty values')
    parser.add_argument("-i", default = None, dest = 'input_file', help="Input file")
    parser.add_argument("-o", default = None, dest = 'output_file', help="Output file")
    parser.add_argument("-t", default = "transcripts.txt", dest = 'transcripts_file', help="Transcripts file")
    args = parser.parse_args()
    main(**vars(args))

if __name__ == '__main__':
    parse()
