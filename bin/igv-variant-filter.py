#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Filters the ANNOVAR annotation .tsv table for usage with IGV Snapshots

INPUT: ANNOVAR annotations merged with original .vcf .tsv table
OUTPUT: Filtered annotations .tsv table
USAGE: igv-variant-filter.py -c HaplotypeCaller -s "sampleID" -i "annotations.tsv" -o "sampleID.tmb.filtered.tsv"


Criteria:

For both matched and unmatched we will apply the following criteria:
1- VAF >5% tumor
2- VAF <2% normal
"""
import csv
import sys
import argparse

from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)
"""
https://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python
"""

frequency_min_tumor = 0.05 # 5%
frequency_min_normal = 0.02 # 5%


def unpaired_filter(row):
    """
    Return True or False if the row passes all the filter criteria
    """
    frequency = float(row['AF'])

    frequency_pass = frequency > frequency_min_tumor

    return(all([ frequency_pass ]))

def LoFreqSomatic(fin, fout):
    reader = csv.DictReader(fin, delimiter = '\t')
    fieldnames = reader.fieldnames
    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = fieldnames)
    writer.writeheader()
    for row in reader:
        if unpaired_filter(row):
            writer.writerow(row)

def MuTect2(fin, fout):
    reader = csv.DictReader(fin, delimiter = '\t')
    fieldnames = reader.fieldnames
    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = fieldnames)
    writer.writeheader()
    for row in reader:
        if unpaired_filter(row):
            writer.writerow(row)



def main(**kwargs):
    """
    Main control function for the script
    """
    input_file = kwargs.pop('input_file', None)
    output_file = kwargs.pop('output_file', None)
    caller = kwargs.pop('caller')

    if input_file:
        fin = open(input_file)
    else:
        fin = sys.stdin

    if output_file:
        fout = open(output_file, "w")
    else:
        fout = sys.stdout

    if caller == "LoFreqSomatic":
        LoFreqSomatic(fin, fout)
        fout.close()
        fin.close()
    elif caller == "MuTect2":
        MuTect2(fin, fout) # TODO: create this function & filter methods for paired calling
        fout.close()
        fin.close()
    else:
        print("ERROR: caller not recognized: {0}".format(caller))
        sys.exit(1)

def parse():
    """
    Parses script args
    """
    parser = argparse.ArgumentParser(description='Filters the ANNOVAR annotation .tsv table for usage with IGV snapshots')
    parser.add_argument("-i", default = None, dest = 'input_file', help="Input file")
    parser.add_argument("-o", default = None, dest = 'output_file', help="Output file")
    parser.add_argument("-c", "--caller", dest = 'caller', help="Variant caller used", required=True)
    args = parser.parse_args()

    main(**vars(args))

if __name__ == '__main__':
    parse()
