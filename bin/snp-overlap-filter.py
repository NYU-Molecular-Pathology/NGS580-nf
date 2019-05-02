#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Filters the ANNOVAR annotation .tsv table for usage with Homozygous SNP Overlap

INPUT: ANNOVAR annotations merged with original .vcf .tsv table
OUTPUT: Filtered annotations .tsv table
USAGE: snp-overlap-filter.py -c HaplotypeCaller -i "annotations.tsv" -o "sampleID.snp-overlap.filtered.tsv"


Filter Criteria:

Coverage >200
Variant Allele Frequency >0.98 (98%)
"""
import csv
import sys
import argparse

from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)
"""
https://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python
"""

coverage_min = 200.0
frequency_min = 0.98 # 98%
# ExonicFunc_refGene_not_allowed_substr = [ 'deletion', 'insertion' ]

def unpaired_filter(row):
    """
    Return True or False if the row passes all the filter criteria
    """
    frequency = float(row['AF'])
    coverage = float(row['DP'])
    ExonicFunc_refGene = row['ExonicFunc.refGene']

    frequency_pass = frequency > frequency_min
    coverage_pass = coverage > coverage_min
    # ExonicFunc_refGene_pass = False
    # for substr in ExonicFunc_refGene_not_allowed_substr:
    #     ExonicFunc_refGene_pass = not substr in ExonicFunc_refGene

    return(all([ coverage_pass, frequency_pass ]))

def LoFreq(fin, fout):
    reader = csv.DictReader(fin, delimiter = '\t')
    fieldnames = reader.fieldnames
    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = fieldnames)
    writer.writeheader()
    for row in reader:
        if unpaired_filter(row):
            writer.writerow(row)

def HaplotypeCaller(fin, fout):
    reader = csv.DictReader(fin, delimiter = '\t')
    fieldnames = reader.fieldnames
    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = fieldnames)
    writer.writeheader()
    for row in reader:
        if unpaired_filter(row):
            writer.writerow(row)

def VarScan2(fin, fout):
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

    if caller == "HaplotypeCaller":
        HaplotypeCaller(fin, fout)
        fout.close()
        fin.close()
    elif caller == "LoFreq":
        LoFreq(fin, fout)
        fout.close()
        fin.close()
    elif caller == "MuTect2":
        MuTect2(fin, fout) # TODO: create this function & filter methods for paired calling
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
    parser = argparse.ArgumentParser(description='Filters the ANNOVAR annotation .tsv table for usage with Homozygous SNP Overlap')
    parser.add_argument("-i", default = None, dest = 'input_file', help="Input file")
    parser.add_argument("-o", default = None, dest = 'output_file', help="Output file")
    parser.add_argument("-c", "--caller", dest = 'caller', help="Variant caller used", required=True)
    args = parser.parse_args()

    main(**vars(args))

if __name__ == '__main__':
    parse()
