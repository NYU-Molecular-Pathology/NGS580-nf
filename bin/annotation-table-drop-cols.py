#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Remove extraneous columns from the annotation table
"""
import csv
import sys
import argparse

drop_cols_all = [
'FREQ',
'NLOD',
'TLOD',
'NORMAL.AD',
'TUMOR.AD',
'FILTER',
'INDEL',
'ID',
'CONSVAR',
'HRUN',
'ExAC_AFR',
'ExAC_AMR',
'ExAC_EAS',
'ExAC_FIN',
'ExAC_NFE',
'ExAC_OTH',
'ExAC_SAS',
'AC',
'AN',
'RD',
'RBQ',
'ABQ',
'QUAL.1',
'END',
'HOMLEN',
'SVLEN',
'SVTYPE',
'SOMATIC',
'QSS',
'MQ',
'SNVSB',
'SomaticEVS',
'NORMAL.AU',
'NORMAL.TU',
'NORMAL.CU',
'NORMAL.GU',
'TUMOR.AU',
'TUMOR.TU',
'TUMOR.CU',
'TUMOR.GU',
'DP4',
'UQ',
'UNIQ',
'QSI',
'NORMAL.TAR',
'NORMAL.TIR',
'NORMAL.TOR',
'TUMOR.TAR',
'TUMOR.TIR',
'TUMOR.TOR',
'CLNDSDBID',
'CLNDSDB',
'CLNACC',
'CLNDBN',
'CLINSIG',
'1000g2015aug_all',
'avsnp150',
'ExAC_ALL',
'snp138',
'cosmic70',
'Chr',
'Start',
'End',
'Ref',
'Alt',
'SB'
]

drop_cols_paired = [
'Sample',
'AD.REF',
'AD.ALT',
'AF.ALT',
'AF.REF'
]

drop_cols_unpaired = [
'NORMAL.AD',
'NORMAL.DP',
'NORMAL.AF',
'TUMOR.AD',
'TUMOR.DP',
'TUMOR.AF',
'TUMOR.AD.REF',
'TUMOR.AD.ALT',
'TUMOR.AD.TOTAL',
'NORMAL.AD.REF',
'NORMAL.AD.ALT',
'NORMAL.AD.TOTAL',
'Tumor',
'Normal'
]

def main(**kwargs):
    """
    Main control function for the script
    """
    input_file = kwargs.pop('input_file', None)
    output_file = kwargs.pop('output_file', None)
    type = kwargs.pop('type', 'all')

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
    old_fieldnames = reader.fieldnames

    # remove the unwanted cols from fieldnames
    new_fieldnames = [ f for f in old_fieldnames if f not in drop_cols_all ]
    if type == 'paired':
        new_fieldnames = [ f for f in new_fieldnames if f not in drop_cols_paired ]
    if type == 'unpaired':
        new_fieldnames = [ f for f in new_fieldnames if f not in drop_cols_unpaired ]

    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = new_fieldnames)
    writer.writeheader()

    for row in reader:
        for key in drop_cols_all:
            row.pop(key, None)
        if type == 'paired':
            for key in drop_cols_paired:
                row.pop(key, None)
        if type == 'unpaired':
            for key in drop_cols_unpaired:
                row.pop(key, None)
        writer.writerow(row)

    fout.close()
    fin.close()


def parse():
    """
    Parses script args
    """
    parser = argparse.ArgumentParser(description='Remove extraneous columns from the annotation table')
    parser.add_argument("-i", default = None, dest = 'input_file', help="Input file")
    parser.add_argument("-o", default = None, dest = 'output_file', help="Output file")
    parser.add_argument("--type", default = 'all', dest = 'type',
        choices=['paired', 'unpaired', 'all'],
        help="Type of table to apply extra column filters on ('paired' or 'unpaired')")
    args = parser.parse_args()
    main(**vars(args))

if __name__ == '__main__':
    parse()
