#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Remove extraneous columns from the annotation table
"""
import csv
import sys
import argparse

drop_cols = [
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
'TUMOR.TOR'
]

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
    old_fieldnames = reader.fieldnames

    # remove the unwanted cols from fieldnames
    new_fieldnames = [ f for f in old_fieldnames if f not in drop_cols ]
    # for fieldname in old_fieldnames:
    #     if fieldname not in drop_cols:
    #         new_fieldnames.append(fieldname)

    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = new_fieldnames)
    writer.writeheader()

    for row in reader:
        for key in drop_cols:
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
    args = parser.parse_args()
    main(**vars(args))

if __name__ == '__main__':
    parse()
