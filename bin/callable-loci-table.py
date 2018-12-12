#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert GATK CallableLoci summary table to a real table format

example input file looks like this:

                         state nBases
                         REF_N 0
                      CALLABLE 1677650
                   NO_COVERAGE 23528
                  LOW_COVERAGE 329335
            EXCESSIVE_COVERAGE 0
          POOR_MAPPING_QUALITY 0
"""
import sys
import csv

inputSummary = sys.argv[1] # "output/Sample1.CallableLoci.summary.txt"
outputTable = sys.argv[2] # "output.tsv"

entries = []

with open(inputSummary) as fin:
    headerLine = fin.readline().strip()
    headers = headerLine.split(' ')
    for line in fin:
        parts = line.strip().split(' ')
        line_dict = {}
        for i, header in enumerate(headers):
            line_dict[header] = parts[i]
        entries.append(line_dict)

with open(outputTable, "w") as fout:
    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = headers)
    writer.writeheader()
    for entry in entries:
        writer.writerow(entry)
