#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Reformats the .tsv vcf table output to recalculate and standardize values for downstream processing.

INPUT: .TSV formatted output from GATK VariantsToTable
OUTPUT: .TSV formatted table with recalculated variant allele frequency in the FREQ column
USAGE: reformat-vcf-table.py -c HaplotypeCaller -s "sampleID" -i "sample_tsv" -o "sampleID.recalc.tsv"
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

def HaplotypeCaller(fin, fout, sampleID):
    """
    Recalculates the variant allele frequency for GATK Haplotype Caller output
    assumes VCFv4.2
    Outputs extra columns using variant values

    Parameters
    ----------
    fin: file connection
        connection to the input file
    fout: file connection
        connection to the output file
    sampleID: str
        identifier for the sample in the input file connection

    """
    # column names for the AD and DP fields in the table
    AD_key = "{0}.AD".format(sampleID)
    DP_key = "{0}.DP".format(sampleID)

    reader = csv.DictReader(fin, delimiter = '\t')
    fieldnames = reader.fieldnames
    fieldnames.append('FREQ')
    fieldnames.append('AD.REF')
    fieldnames.append('AD.ALT')
    fieldnames.append('AF.ALT')
    fieldnames.append('AF.REF')
    fieldnames.append('DP')
    fieldnames.append('AF')
    fieldnames_out = [ item for item in fieldnames ]
    fieldnames_out.remove(AD_key)
    fieldnames_out.remove(DP_key)
    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = fieldnames_out)
    writer.writeheader()
    for row in reader:
        ref_AD = float(row[AD_key].split(',')[0])
        alt_AD = float(row[AD_key].split(',')[1])
        depth = float(row[DP_key])
        ref_AF = ref_AD / depth
        alt_AF = alt_AD / depth
        row['FREQ'] = alt_AF
        row['AF'] = alt_AF
        row['AD.REF'] = ref_AD
        row['AD.ALT'] = alt_AD
        row['AF.REF'] = ref_AF
        row['AF.ALT'] = alt_AF
        row['DP'] = depth
        drop_vals = []
        drop_vals.append(row.pop(AD_key))
        drop_vals.append(row.pop(DP_key))
        writer.writerow(row)

def LoFreq(fin, fout):
    """
    LoFreq does not need recalulating; output a new column in the file with 'FREQ' for consistency

    Parameters
    ----------
    fin: file connection
        connection to the input file
    fout: file connection
        connection to the output file

    """
    # allele frequency column
    AF_key = "AF"
    reader = csv.DictReader(fin, delimiter = '\t')
    fieldnames = reader.fieldnames
    fieldnames.append('FREQ')
    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = fieldnames)
    writer.writeheader()
    for row in reader:
        AF_val = row[AF_key]
        row['FREQ'] = AF_val
        writer.writerow(row)

def MuTect2(fin, fout):
    """
    Outputs the TUMOR.AF as a new column called 'FREQ', and outputs the TLOD value as QUAL
    Adds extra allelic depth columns

    Parameters
    ----------
    fin: file connection
        connection to the input file
    fout: file connection
        connection to the output file

    """
    reader = csv.DictReader(fin, delimiter = '\t')
    # get old headers
    fieldnames = reader.fieldnames
    # append new headers for the columns to be created
    fieldnames.append('FREQ')
    fieldnames.append('QUAL')
    fieldnames.append('TUMOR.AD.REF')
    fieldnames.append('TUMOR.AD.ALT')
    fieldnames.append('TUMOR.AD.TOTAL')
    fieldnames.append('NORMAL.AD.REF')
    fieldnames.append('NORMAL.AD.ALT')
    fieldnames.append('NORMAL.AD.TOTAL')
    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = fieldnames)
    writer.writeheader()
    for row in reader:
        # set the FREQ to the tumor AF
        tumor_AF_value = row["TUMOR.AF"]
        row['FREQ'] = row['TUMOR.AF']
        # change QUAL to the TLOD
        row['QUAL'] = row['TLOD']
        # split the tumor AD values for ref and alt
        row['TUMOR.AD.REF'] = row['TUMOR.AD'].split(',')[0]
        row['TUMOR.AD.ALT'] = row['TUMOR.AD'].split(',')[1]
        row['NORMAL.AD.REF'] = row['NORMAL.AD'].split(',')[0]
        row['NORMAL.AD.ALT'] = row['NORMAL.AD'].split(',')[1]
        # add up the total allelic depths
        row['TUMOR.AD.TOTAL'] = int(row['TUMOR.AD.REF']) + int(row['TUMOR.AD.ALT'])
        row['NORMAL.AD.TOTAL'] = int(row['NORMAL.AD.REF']) + int(row['NORMAL.AD.ALT'])
        writer.writerow(row)

def VarScan2(fin, fout, sampleID):
    """
    Reformat the contents of lines VarScan2 output
    VarScan2 outputs the variant allele frequency ('FREQ') like this: "100%", "99.1%", "53.42%", etc
    Need to convert this to decimal float

    Many of the table headers in VarScan2 output have the SampleID prepended to them

    - need to rename the FREQ columns from "Sample1.FREQ" to "FREQ"
    - change FREQ value from character string to float; 99.91% -> 0.99
    - need to add AF column, set to value of FREQ

    NOTE: Make sure this corresponds to the columns output in the vcf_to_tsv pipeline step by GATK VariantsToTable
    """
    reader = csv.DictReader(fin, delimiter = '\t')

    # get old headers
    # ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'Sample1.AD', 'Sample1.RD', 'Sample1.FREQ', 'Sample1.RBQ', 'Sample1.ABQ']
    old_fieldnames = reader.fieldnames

    # get the headers that dont have sample ID embedded in them
    # ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']
    new_fieldnames = [ n for n in old_fieldnames if sampleID not in n ]
    # add the missing field names
    new_fieldnames.append('DP')
    new_fieldnames.append('AD')
    new_fieldnames.append('RD')
    new_fieldnames.append('FREQ')
    new_fieldnames.append('AF')
    new_fieldnames.append('RBQ')
    new_fieldnames.append('ABQ')
    new_fieldnames.append('QUAL.REF')
    new_fieldnames.append('QUAL.ALT')
    new_fieldnames.append('AD.ALT')
    new_fieldnames.append('AD.REF')

    old_DP = "{0}.DP".format(sampleID)
    old_AD = "{0}.AD".format(sampleID)
    old_RD = "{0}.RD".format(sampleID)
    old_FREQ = "{0}.FREQ".format(sampleID)
    old_RBQ = "{0}.RBQ".format(sampleID)
    old_ABQ = "{0}.ABQ".format(sampleID)

    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = new_fieldnames)
    writer.writeheader()
    for row in reader:
        # convert the old row columns to the new header keys
        # fix the FREQ value
        # strip the percent sign
        row['FREQ'] = re.sub('[%]', '', row[old_FREQ])
        # divide by 100
        row['FREQ'] = float(row['FREQ']) / 100.0
        # truncate to two decimal places
        row['FREQ'] = '{:0.2f}'.format(row['FREQ'])
        row['AF'] = row['FREQ']

        # fill in the missing required columns
        row['DP'] = row[old_DP]
        row['AD'] = row[old_AD]
        row['RD'] = row[old_RD]
        row['RBQ'] = row[old_RBQ]
        row['ABQ'] = row[old_ABQ]

        row['QUAL'] = row['ABQ']
        row['QUAL.ALT'] = row['ABQ']
        row['QUAL.REF'] = row['RBQ']

        row['AD.ALT'] = row['AD']
        row['AD.REF'] = row['RD']

        # get rid of the columns with old bad headers
        row.pop(old_DP)
        row.pop(old_AD)
        row.pop(old_RD)
        row.pop(old_RBQ)
        row.pop(old_ABQ)
        row.pop(old_FREQ)

        writer.writerow(row)


def main(**kwargs):
    """
    Main control function for the script
    """
    input_file = kwargs.pop('input_file', None)
    output_file = kwargs.pop('output_file', None)
    caller = kwargs.pop('caller')
    sampleID = kwargs.pop('sampleID')

    if input_file:
        fin = open(input_file)
    else:
        fin = sys.stdin

    if output_file:
        fout = open(output_file, "w")
    else:
        fout = sys.stdout

    if caller == "HaplotypeCaller":
        HaplotypeCaller(fin, fout, sampleID)
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
        VarScan2(fin, fout, sampleID)
        fout.close()
        fin.close()
    else:
        print("ERROR: caller not recognized: {0}".format(caller))
        sys.exit(1)




def parse():
    """
    Parses script args
    """
    parser = argparse.ArgumentParser(description='Append a column of text to a file')
    parser.add_argument("-i", default = None, dest = 'input_file', help="Input file")
    parser.add_argument("-o", default = None, dest = 'output_file', help="Output file")

    parser.add_argument("-c", "--caller", dest = 'caller', help="Variant caller used", required=True)
    parser.add_argument("-s", "--sampleID", dest = 'sampleID', help="Sample ID", required=True)
    args = parser.parse_args()

    main(**vars(args))



if __name__ == '__main__':
    parse()
