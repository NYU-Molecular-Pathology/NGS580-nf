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

def filter_row(row, seracare_selected):
    """
    Return True or False if the row passes all the filter criteria
    """
    allowed = False
    key = "{0}{1}{2}{3}".format(row['CHROM'], row['POS'], row['REF'], row['ALT'])
    if key in seracare_selected:
        allowed = True
        row['SeraCare.Gene'] = seracare_selected[key]['Gene']
        row['SeraCare.Coding'] = seracare_selected[key]['Coding']
        row['SeraCare.COSMIC'] = seracare_selected[key]['COSMIC']
        row['SeraCare.AAChange'] = seracare_selected[key]['AAChange']
        row['SeraCare.AF'] = seracare_selected[key]['AF']

    return(allowed, row)

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
    seracare_selected_list = []
    with open(seracare_tsv) as f:
        reader = csv.DictReader(f, delimiter = '\t')
        for row in reader:
            seracare_selected_list.append(row)
    # convert to a dict
    seracare_selected = {}
    for item in seracare_selected_list:
        key = "{0}{1}{2}{3}".format(item['CHROM'], item['POS'], item['REF'], item['ALT'])
        seracare_selected[key] = item

    # start loading the sample annotations
    reader = csv.DictReader(fin, delimiter = '\t')
    # need to update the output fieldnames to include the SeraCare Selected columns
    fieldnames = reader.fieldnames
    fieldnames.append('SeraCare.Gene')
    fieldnames.append('SeraCare.Coding')
    fieldnames.append('SeraCare.AAChange')
    fieldnames.append('SeraCare.COSMIC')
    fieldnames.append('SeraCare.AF')
    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = fieldnames)
    writer.writeheader()
    for row in reader:
        # update row columns with default values
        row['SeraCare.Gene'] = '.'
        row['SeraCare.Coding'] = '.'
        row['SeraCare.COSMIC'] = '.'
        row['SeraCare.AAChange'] = '.'
        row['SeraCare.AF'] = '.'
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
