#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Splits the annotations.tsv file into separate files per variant caller
"""
import sys
import csv

def main(input_file):
    """
    """
    # get all the variant callers
    all_variant_callers = []
    with open(input_file) as fin:
        first_line = fin.readline()
        variant_caller_column = first_line.split('\t').index("VariantCaller")
        for line in fin:
            variant_caller = line.split('\t')[variant_caller_column].strip()
            all_variant_callers.append(variant_caller)
    all_variant_callers = list(set(all_variant_callers))
    
    # set the output files
    variant_caller_outputs = {caller: "annotations." + caller + ".tsv" for caller in all_variant_callers }
    
    # write header to each output file
    for output_file in variant_caller_outputs.values():
        with open(output_file, "w") as fout:
            fout.write(first_line)
    
    # write the remaining lines to all files
    with open(input_file) as fin:
        first_line = fin.readline()
        for line in fin:
            variant_caller = line.split('\t')[variant_caller_column].strip()
            output_file = variant_caller_outputs[variant_caller]
            with open(output_file, "a") as fout:
                fout.write(line)

def parse():
    """
    Parse script args if called as script
    """
    input_file = sys.argv[1]
    main(input_file)

if __name__ == '__main__':
    parse()
