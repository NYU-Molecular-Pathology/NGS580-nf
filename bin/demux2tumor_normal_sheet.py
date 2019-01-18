#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert a demultiplexing samplesheet to a samples.tumor.normal.csv sheet

Usage:
bin/demux2tumor_normal_sheet.py SampleSheet.csv samples.tumor.normal.csv
"""
import sys
import csv
from util import samplesheet

samplesheet_file = sys.argv[1]
output_file = sys.argv[2]
tumor_colname = 'Tumor'
normal_colname = 'Normal'

# load demultiplexing samplesheet
sheet_obj = samplesheet.IEMFile(path = samplesheet_file)
sheet_obj.isValid(_raise = True)

tumor_normal_samples = []

# get all the paired samples from the samplesheet
for sample in sheet_obj.data['Data']['Samples']:
    if 'Paired_Normal' in sample:
        paired_normal = sample['Paired_Normal'].strip()
        if len(sample['Paired_Normal'].strip()) > 0:
            sample_dict = {
            tumor_colname: sample['Sample_ID'],
            normal_colname: sample['Paired_Normal']
            }
            tumor_normal_samples.append(sample_dict)

# write new sheet with just the paired samples
with open(output_file, 'w') as f:
    writer = csv.DictWriter(f, fieldnames = [tumor_colname, normal_colname])
    writer.writeheader()
    for item in tumor_normal_samples:
        writer.writerow(item)
