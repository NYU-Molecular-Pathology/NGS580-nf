#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calcalate the tumor mutation burden in variants/Megabase
"""
import os
import csv
import argparse

"""
Additional criteria:

For both matched and unmatched we will apply the following criteria:
1- VAF >5%
2- coverage > 200X.
3- include SNV only including synonymous and Non Synonymous.
4- exclude Cosmic HS (not apply)
5- include only exonic.

For Unmatched. Add condition below;
6- Frequency of <0.4% in exAc.
"""

NA_strs = ['.'] # strings used as "NA" values in the table
Func_refGene_allowed = ['exonic'] # , 'splicing' 'exonic;splicing', , 'UTR5'
coverage_min = 200.0 # should correspond to GATK CallableLoci depth cutoff
frequency_min = 0.05
ExAC_max = 0.004
ExonicFunc_refGene_allowed = ['synonymous SNV','nonsynonymous SNV']

def filter_rules(row):
    """
    Return True or False if the row passes all the filter criteria
    """
    frequency = float(row['AF'])
    coverage = float(row['DP'])
    #COSMIC = row['cosmic70']
    Func_refGene = row['Func.refGene']
    ExonicFunc_refGene = row['ExonicFunc.refGene']
    ExAC_value = row['ExAC_ALL']

    frequency_pass = frequency > frequency_min
    coverage_pass = coverage >= coverage_min
    #not_in_COSMIC = COSMIC in NA_strs
    not_in_COSMIC = True
    in_Func_refGene_allowed = Func_refGene in Func_refGene_allowed
    in_Func_ExonicFunc_refGene = ExonicFunc_refGene in ExonicFunc_refGene_allowed
    in_ExAC_allowed = ExAC_value in NA_strs or float(ExAC_value) <= ExAC_max
    in_ExAC_allowed = True
    return(all([in_Func_refGene_allowed, not_in_COSMIC, coverage_pass,
                frequency_pass, in_ExAC_allowed, in_Func_ExonicFunc_refGene]))

def get_options():

    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--loci", type=str, required=True,
                        help="callable loci from alignments")
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Path to input variant annotation including file name")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Path to output TMB file including file name")
    return parser.parse_args()

def main():
    args = get_options()
    print(args)
    loci_dict = dict()
    tmb_dict = dict()
    with open (args.loci) as fl:
        for line in fl:
            items = line.strip().split()
            loci_dict[items[1]] = items[0]
    with open(args.input) as fin, open(args.output, 'w') as fout:
        fout.write("SampleID\tVariantCaller\tnBases\tnVariants\tTMB\n")
        reader = csv.DictReader(fin, delimiter = '\t')
        print(reader.fieldnames)
        for row in reader:
            if filter_rules(row):
                caller = row['VariantCaller']
                if caller not in tmb_dict.keys():
                    tmb_dict[caller] = {}
                sample = row['Tumor']
                if sample not in tmb_dict[caller].keys():
                    tmb_dict[caller][sample] = {"variants":0}
                tmb_dict[caller][sample]['variants'] += 1

        print (tmb_dict)
        for caller in tmb_dict.keys():
            for sample in tmb_dict[caller].keys():
                val = tmb_dict[caller][sample]
                tmb = round(float(val['variants'])/float(loci_dict[sample])*1000000,2)
                fout.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(sample,caller,str(loci_dict[sample]),
                                                    str(val['variants']),str(tmb), caller))

if __name__ == "__main__":
    main()
