#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Filter rows from the final merged annotation table, to apply filter
criteria that could not be easily applied earlier in the pipeline
"""
import csv
import sys
import argparse

lofreq_min_frequency = 0.02
lofreq_max_frequency = 0.99
depth_min = 100
somatic_frequency_min_tumor = 0.05 # 5% tumor AF
somatic_frequency_min_normal = 0.02 # 2% normal AF
Func_refGene_exclude = [
'ncRNA_intronic',
'intronic',
'UTR3'
]

def HaplotypeCaller(row):
    return(True)
    funcRefGene = row['Func.refGene']
    passFuncRefGene = True
    if funcRefGene in Func_refGene_exclude:
        passFuncRefGene = False

    if passFuncRefGene:
        return(True)
    else:
        return(False)

def LoFreq(row):
    return(True)
    # depth = float(row['DP'])
    # alleleFrequency = float(row['AF'])
    #
    # passMinAF = True
    # passMaxAF = True
    # if alleleFrequency < lofreq_min_frequency:
    #     passMinAF = False
    # if alleleFrequency > lofreq_max_frequency:
    #     passMaxAF = False
    #
    # depthPass = True
    # if depth < depth_min:
    #     depthPass = False

    funcRefGene = row['Func.refGene']
    passFuncRefGene = True
    if funcRefGene in Func_refGene_exclude:
        passFuncRefGene = False

    if passFuncRefGene: # passMinAF and passMaxAF and depthPass and
        return(True)
    else:
        return(False)

def MuTect2(row):
    # tumorAlleleFrequency = float(row['AF'])
    # if tumorAlleleFrequency > somatic_frequency_min_tumor:
    #     return(True)
    # else:
    #     return(False)
    return(True)

def VarScan2(row):
    return(True)
    # do not include VarScan2 output right now
    return(False)

def StrelkaSomaticIndel(row):
    return(True)
    # AF not encoded in .vcf
    tumorAlleleFrequency = float(row['AF'])
    passAF = True
    if tumorAlleleFrequency < somatic_frequency_min_tumor:
        passAF = False

    funcRefGene = row['Func.refGene']
    passFuncRefGene = True
    if funcRefGene in Func_refGene_exclude:
        passFuncRefGene = False

    if passFuncRefGene and passAF:
        return(True)
    else:
        return(False)

def StrelkaSomaticSNV(row):
    return(True)
    tumorAlleleFrequency = float(row['AF'])
    passAF = True
    if tumorAlleleFrequency < somatic_frequency_min_tumor:
        passAF = False

    funcRefGene = row['Func.refGene']
    passFuncRefGene = True
    if funcRefGene in Func_refGene_exclude:
        passFuncRefGene = False

    if passFuncRefGene and passAF:
        return(True)
    else:
        return(False)

def Pindel(row):
    return(True)
    funcRefGene = row['Func.refGene']
    passFuncRefGene = True
    if funcRefGene in Func_refGene_exclude:
        passFuncRefGene = False

    # AF not encoded in vcf
    tumorAlleleFrequency = float(row['AF'])
    passAF = True
    if tumorAlleleFrequency < somatic_frequency_min_tumor:
        passAF = False

    # DP not encoded in vcf
    depth = float(row['DP'])
    depthPass = True
    if depth < depth_min:
        depthPass = False

    if passAF and passFuncRefGene and depthPass:
        return(True)
    else:
        return(False)

def LoFreqSomatic(row):
    return(True)
    # tumorAlleleFrequency = float(row['AF'])
    # passMinAF = True
    # passMaxAF = True
    # if tumorAlleleFrequency < lofreq_min_frequency:
    #     passMinAF = False
    # if tumorAlleleFrequency > lofreq_max_frequency:
    #     passMaxAF = False

    funcRefGene = row['Func.refGene']
    passFuncRefGene = True
    if funcRefGene in Func_refGene_exclude:
        passFuncRefGene = False

    if passFuncRefGene: # and passMinAF and passMaxAF
        return(True)
    else:
        return(False)


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
    fieldnames = reader.fieldnames
    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = fieldnames)
    writer.writeheader()

    for row in reader:
        variantCaller = row['VariantCaller']
        variantCallerType = row['VariantCallerType']

        if variantCaller == "HaplotypeCaller":
            if HaplotypeCaller(row):
                writer.writerow(row)

        elif variantCaller == "LoFreq":
            if LoFreq(row):
                writer.writerow(row)

        elif variantCaller == "MuTect2":
            if MuTect2(row):
                writer.writerow(row)

        elif variantCaller == "VarScan2":
            if VarScan2(row):
                writer.writerow(row)

        elif variantCaller == "Strelka":
            if variantCallerType == "snvs":
                if StrelkaSomaticIndel(row):
                    writer.writerow(row)
            if variantCallerType == "indel":
                if StrelkaSomaticSNV(row):
                    writer.writerow(row)

        elif variantCaller == "Pindel":
            if Pindel(row):
                writer.writerow(row)

        elif variantCaller == "LoFreqSomatic":
            if LoFreqSomatic(row):
                writer.writerow(row)

        else:
            print("ERROR: variantCaller not recognized: {0}".format(variantCaller))
            sys.exit(1)

    fout.close()
    fin.close()


def parse():
    """
    Parses script args
    """
    parser = argparse.ArgumentParser(description='Filter rows from the annotation table')
    parser.add_argument("-i", default = None, dest = 'input_file', help="Input file")
    parser.add_argument("-o", default = None, dest = 'output_file', help="Output file")
    args = parser.parse_args()
    main(**vars(args))

if __name__ == '__main__':
    parse()
