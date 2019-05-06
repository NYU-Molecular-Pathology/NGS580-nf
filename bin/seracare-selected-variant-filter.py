#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Filters the ANNOVAR annotation .tsv table for variants present in the SeraCare Selected Variants .tsv file
"""
import csv
import sys
import argparse

from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)
"""
https://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python
"""

# NA_strs = ['.'] # strings used as "NA" values in the table
# Func_refGene_allowed = ['exonic'] # , 'splicing' 'exonic;splicing', , 'UTR5'
# coverage_min = 500.0 # should correspond to GATK CallableLoci depth cutoff
# frequency_min = 0.05 # 5%
# ExAC_allowed = ['.', '0'] # only allow NA or 0 values

def filter_row(row, seracare_selected):
    """
    Return True or False if the row passes all the filter criteria
    """
    # update row columns with default values
    row['SeraCare.Gene'] = '.'
    row['SeraCare.COSMIC'] = '.'
    row['SeraCare.AminoAcid'] = '.'
    row['SeraCare.Type'] = '.'
    row['SeraCare.TargetAF'] = '.'

    # COSMIC id annotation string from ANNOVAR
    row_COSMIC = row['cosmic70']

    # false by default
    in_allowed_COSMIC = False

    # check if the COSMIC id in the sample annotation is included in the SeraCare Selected COSMIC identifiers
    for selected in seracare_selected:
        if selected['COSMIC'] in row_COSMIC:
            in_allowed_COSMIC = True
            row['SeraCare.Gene'] = selected['Gene']
            row['SeraCare.COSMIC'] = selected['COSMIC']
            row['SeraCare.AminoAcid'] = selected['AminoAcid']
            row['SeraCare.Type'] = selected['Type']
            row['SeraCare.TargetAF'] = selected['TargetAF']

    return(in_allowed_COSMIC, row)

def main(**kwargs):
    """
    Main control function for the script
    """
    input_file = kwargs.pop('input_file', None)
    output_file = kwargs.pop('output_file', None)
    seracare_tsv = kwargs.pop('seracare_tsv')

    if input_file:
        fin = open(input_file)
    else:
        fin = sys.stdin

    if output_file:
        fout = open(output_file, "w")
    else:
        fout = sys.stdout

    # load the SeraCare Selected Variants from file
    seracare_selected = []
    with open(seracare_tsv) as f:
        reader = csv.DictReader(f, delimiter = '\t')
        for row in reader:
            seracare_selected.append(row)

    # start loading the sample annotations
    reader = csv.DictReader(fin, delimiter = '\t')
    # need to update the output fieldnames to include the SeraCare Selected columns
    fieldnames = reader.fieldnames
    fieldnames.append('SeraCare.Gene')
    fieldnames.append('SeraCare.COSMIC')
    fieldnames.append('SeraCare.AminoAcid')
    fieldnames.append('SeraCare.Type')
    fieldnames.append('SeraCare.TargetAF')
    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = fieldnames)
    writer.writeheader()
    for row in reader:
        allowed, new_row = filter_row(row = {k:v for k,v in row.items()}, seracare_selected = seracare_selected)
        if allowed:
            writer.writerow(new_row)

def parse():
    """
    Parses script args
    """
    parser = argparse.ArgumentParser(description='Filters the ANNOVAR annotation .tsv table for usage with Tumor Mutation Burden analysis')
    parser.add_argument("-i", default = None, dest = 'input_file', help="Input file")
    parser.add_argument("-o", default = None, dest = 'output_file', help="Output file")
    parser.add_argument("-s", "--seracare", required = True, dest = 'seracare_tsv', help="SeraCare Selected Variants TSV")
    args = parser.parse_args()

    main(**vars(args))

if __name__ == '__main__':
    parse()
