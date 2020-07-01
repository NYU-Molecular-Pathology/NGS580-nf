#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Splits the annotations.tsv file into separate files
based on variant callers
and whether the variant caller is designated at being used for
"paired" or "unpaired" variant calling (tumor normal pairs)
"""
import sys

def main(input_file, output_paired, output_unpaired):
    """
    """
    paired_callers = [
    'LoFreqSomatic',
    'MuTect2',
    'Strelka'
    ]
    fin = open(input_file)
    paired_out = open(output_paired, "w")
    unpaired_out = open(output_unpaired, "w")

    # get and print the header
    first_line = fin.readline()
    unpaired_out.write(first_line)
    paired_out.write(first_line)

    # get index of the column that has the VariantCaller
    variant_caller_column_index = first_line.split('\t').index("VariantCaller")

    # iterate over remaining lines and write lines to the different outputs
    for line in fin:
        variant_caller = line.split('\t')[variant_caller_column_index].strip()
        if variant_caller in paired_callers:
            paired_out.write(line)
        else:
            unpaired_out.write(line)

    unpaired_out.close()
    paired_out.close()
    fin.close()


def parse():
    """
    Parse script args if called as script
    """
    input_file = sys.argv[1]
    output_paired = sys.argv[2]
    output_unpaired = sys.argv[3]
    main(input_file, output_paired, output_unpaired)

if __name__ == '__main__':
    parse()
