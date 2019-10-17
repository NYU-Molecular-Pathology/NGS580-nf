#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Split rows in the annotation table that have multiple transcripts listed
So that each transcript is output on its own row

Example:

MNX1:NM_001165255:exon2:c.179_180insGAG:p.F60delinsLS,MNX1:NM_005515:exon2:c.815_816insGAG:p.F272delinsLS

MNX1:NM_001165255:exon2:c.179_180insGAG:p.F60delinsLS
MNX1:NM_005515:exon2:c.815_816insGAG:p.F272delinsLS
"""
import csv
import sys
import argparse

transcript_colname = "AAChange.refGene"
NA_char = "."
invalid_chars = ["UNKNOWN", NA_char]

# MNX1:NM_001165255:exon2:c.179_180insGAG:p.F60delinsLS
transcript_order = {
0 : "Gene",
1 : 'Transcript',
2 : 'Exon',
3 : 'CodingChange',
4 : 'AAChange'
}

def main(**kwargs):
    """
    Main control function for the script
    """
    input_file = kwargs.pop('input_file', None)
    output_file = kwargs.pop('output_file', None)

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
    fieldnames.remove(transcript_colname)
    fieldnames.append("Gene")
    fieldnames.append("Transcript")
    fieldnames.append("Exon")
    fieldnames.append("CodingChange")
    fieldnames.append("AAChange")
    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = fieldnames)
    writer.writeheader()

    for row in reader:
        # initialize with default values because not all items may be present
        transcript_output = {
        'Gene' : NA_char,
        'Transcript': NA_char,
        'Exon' : NA_char,
        'CodingChange': NA_char,
        'AAChange': NA_char
        }

        if row[transcript_colname] in invalid_chars:
            row.pop(transcript_colname)
            row.update(transcript_output)
            writer.writerow(row)
        else:
            # get all the transcripts in the entry
            line_parts = row[transcript_colname].split(',')
            for line_part in line_parts:
                new_row = {k : v for k, v in row.items()}
                # split the transcript entry apart
                transcript_parts = line_part.split(':')
                for i, part in enumerate(transcript_parts):
                    # update the output with the new value
                    transcript_output[transcript_order[i]] = part
                new_row.pop(transcript_colname)
                new_row.update(transcript_output)
                writer.writerow(new_row)

    fout.close()
    fin.close()


def parse():
    """
    Parses script args
    """
    parser = argparse.ArgumentParser(description='Split rows in the annotation table so that entries with multiple transcripts get output each on a new row')
    parser.add_argument("-i", default = None, dest = 'input_file', help="Input file")
    parser.add_argument("-o", default = None, dest = 'output_file', help="Output file")
    args = parser.parse_args()
    main(**vars(args))

if __name__ == '__main__':
    parse()
