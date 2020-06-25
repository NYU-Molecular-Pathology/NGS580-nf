#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Filters the ANNOVAR annotation .tsv table for usage with Tumor Mutation Burden analysis

INPUT: ANNOVAR annotations merged with original .vcf .tsv table
OUTPUT: Filtered annotations .tsv table
USAGE: tmb-variant-filter.py -c HaplotypeCaller -s "sampleID" -i "annotations.tsv" -o "sampleID.tmb.filtered.tsv"



Filtering methods adapted from:
https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-017-0424-2
https://doi.org/10.1186/s13073-017-0424-2
Analysis of 100,000 human cancer genomes reveals the landscape of tumor mutational burden. Chalmers ZR et. al.
Non-coding alterations were not counted.
Alterations listed as known somatic alterations in COSMIC and truncations in tumor suppressor genes were not counted, since our assay genes are biased toward genes with functional mutations in cancer [63].
Alterations predicted to be germline by the somatic-germline-zygosity algorithm were not counted [64].
Alterations that were recurrently predicted to be germline in our cohort of clinical specimens were not counted.
Known germline alterations in dbSNP were not counted.
Germline alterations occurring with two or more counts in the ExAC database were not counted [65].


Additional criteria:

For both matched and unmatched we will apply the following criteria:
1- VAF >5%
2- coverage > 200X.
3- include SNV only including synonymous ans Non Synonymous.
4- exclude Cosmic HS
5- include only exonic and 5’ UTR.

For Unmatched. Add condition below;
6- Frequency of <1% in exAc.

"""
import csv
import sys
import argparse

from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)
"""
https://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python
"""

NA_strs = ['.'] # strings used as "NA" values in the table
Func_refGene_allowed = ['exonic'] # , 'splicing' 'exonic;splicing', , 'UTR5'
coverage_min = 500.0 # should correspond to GATK CallableLoci depth cutoff
frequency_min = 0.05 # 5%
#ExAC_allowed = ['.', '0'] # only allow NA or 0 values
ExAC_max = 0.004
ExonicFunc_refGene_allowed = ['synonymous SNV','nonsynonymous SNV']

def unpaired_filter(row):
    """
    Return True or False if the row passes all the filter criteria
    """
    frequency = float(row['AF'])
    coverage = float(row['DP'])
    COSMIC = row['cosmic70']
    Func_refGene = row['Func.refGene']
    ExonicFunc_refGene = row['ExonicFunc.refGene']
    ExAC_value = row['ExAC_ALL']

    frequency_pass = frequency > frequency_min
    coverage_pass = coverage > coverage_min
    not_in_COSMIC = COSMIC in NA_strs
    in_Func_refGene_allowed = (Func_refGene in Func_refGene_allowed and ExonicFunc_refGene in ExonicFunc_refGene_allowed)
    in_ExAC_allowed = ExAC_value <= ExAC_max
    return(all([in_Func_refGene_allowed, not_in_COSMIC, coverage_pass, frequency_pass, in_ExAC_allowed]))
    # print(frequency, frequency_pass, coverage, coverage_pass, COSMIC, not_in_COSMIC, Func_refGene, in_Func_refGene_allowed)

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
    parser = argparse.ArgumentParser(description='Filters the ANNOVAR annotation .tsv table for usage with Tumor Mutation Burden analysis')
    parser.add_argument("-i", default = None, dest = 'input_file', help="Input file")
    parser.add_argument("-o", default = None, dest = 'output_file', help="Output file")
    parser.add_argument("-c", "--caller", dest = 'caller', help="Variant caller used", required=True)
    args = parser.parse_args()

    main(**vars(args))

if __name__ == '__main__':
    parse()
