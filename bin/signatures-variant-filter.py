#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Filters the variant table for usage with genomic signatures

INPUT: .TSV formatted output from GATK VariantsToTable
OUTPUT: .TSV formatted table with recalculated variant allele frequency in the FREQ column
USAGE: signatures-variant-filter.py -c HaplotypeCaller -s "sampleID" -i "sample_tsv" -o "sampleID.recalc.tsv"
"""
import csv
import sys
import argparse
import re

from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)
"""
https://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python
"""
def HaplotypeCaller(fin, fout):
    """
    Filter a LoFreq vcf .tsv table for usage with genomic signatures

    Criteria for variant inclusion in output:
    VAF >5%
    coverage > 200X.
    """
    reader = csv.DictReader(fin, delimiter = '\t')
    fieldnames = reader.fieldnames
    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = fieldnames)
    writer.writeheader()
    for row in reader:
        depth = float(row['DP'])
        frequency = float(row['AF'])
        if depth > 200 and frequency > 0.05:
            writer.writerow(row)

def LoFreq(fin, fout):
    """
    Filter a LoFreq vcf .tsv table for usage with genomic signatures

    Criteria for variant inclusion in output:
    VAF >5%
    coverage > 200X.
    """
    reader = csv.DictReader(fin, delimiter = '\t')
    fieldnames = reader.fieldnames
    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = fieldnames)
    writer.writeheader()
    for row in reader:
        depth = int(row['DP'])
        frequency = float(row['AF'])
        if depth > 200 and frequency > 0.05:
            writer.writerow(row)


def VarScan2(fin, fout):
    """
    Filter a VarScan2 vcf .tsv table for usage with genomic signatures

    Criteria for variant inclusion in output:
    VAF >5%
    coverage > 200X.
    """
    reader = csv.DictReader(fin, delimiter = '\t')
    fieldnames = reader.fieldnames
    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = fieldnames)
    writer.writeheader()
    for row in reader:
        depth = int(row['DP'])
        frequency = float(row['FREQ'])
        if depth > 200 and frequency > 0.05:
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

    if caller == "HaplotypeCaller":
        HaplotypeCaller(fin, fout)
        fout.close()
        fin.close()
    elif caller == "LoFreq":
        LoFreq(fin, fout)
        fout.close()
        fin.close()
    elif caller == "MuTect2":
        MuTect2(fin, fout)
        fout.close()
        fin.close()
    elif caller == "VarScan2":
        VarScan2(fin, fout)
        fout.close()
        fin.close()
    else:
        print("ERROR: caller not recognized: {0}".format(caller))
        sys.exit(1)




def parse():
    """
    Parses script args
    """
    parser = argparse.ArgumentParser(description='Filters the variant table for usage with genomic signatures')
    parser.add_argument("-i", default = None, dest = 'input_file', help="Input file")
    parser.add_argument("-o", default = None, dest = 'output_file', help="Output file")
    parser.add_argument("-c", "--caller", dest = 'caller', help="Variant caller used", required=True)
    args = parser.parse_args()

    main(**vars(args))



if __name__ == '__main__':
    parse()
