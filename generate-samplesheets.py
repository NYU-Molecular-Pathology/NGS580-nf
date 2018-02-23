#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script will search the provided directories for fastq.gz files matching naming criteria
and generate sample sheets to use for analysis, in .tsv and .json format.
Multiple directories can be supplied, and only the unique files will be used.

Script overview:
- search for all .fastq.gz files

Usage
-----
Example usage:

    ./gather-fastqs2.py example-data example-data example-data

Output
------
Files output by this script:

    - ``samples.fastq.tsv``: one line per R1 R2 fastq file pair
    - ``samples.fastq.json``: one entry per R1 R2 fastq file pair
    - ``samples.analysis.tsv``: one line per sample with R1 R2 file pairs and metadata
    - ``samples.analysis.json``: one entry per sample with R1 R2 file pairs and metadata

Notes
------
Old samplesheets with the same name will be overwritten

Only file names as described below are supported:

    https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm
    Naming Convention

    FASTQ files are named with the sample name and the sample number, which is a numeric assignment based on the order that the sample is listed in the sample sheet. For example: Data\Intensities\BaseCalls\SampleName_S1_L001_R1_001.fastq.gz
    	▶ 	SampleName—The sample name provided in the sample sheet. If a sample name is not provided, the file name includes the sample ID, which is a required field in the sample sheet and must be unique.
    	▶ 	S1—The sample number based on the order that samples are listed in the sample sheet starting with 1. In this example, S1 indicates that this sample is the first sample listed in the sample sheet.

    Note

    Reads that cannot be assigned to any sample are written to a FASTQ file for sample number 0, and excluded from downstream analysis.
    	▶ 	L001—The lane number.
    	▶ 	R1—The read. In this example, R1 means Read 1. For a paired-end run, there is at least one file with R2 in the file name for Read 2. When generated, index reads are I1 or I2.
    	▶ 	001—The last segment is always 001.
"""
import os
import sys
from bin import find
from bin import collapse
import re
import csv
import json
import argparse

# ~~~~~ FUNCTIONS ~~~~~ #
def main(search_dirs, output_prefix = None, NA_value = "NA", tumor_colname = 'Tumor', normal_colname = 'Normal', r1_colname = 'R1', r2_colname = 'R2', sample_colname = 'Sample'):
    """
    Main control function for the script

    Parameters
    ----------
    search_dirs: list
        list of paths to directories to use
    """
    # validate inputs
    if len(search_dirs) < 1:
        print("ERROR: no directories provided")
        sys.exit(1)
    for search_dir in search_dirs:
        if not os.path.isdir(search_dir):
            print("ERROR: '{0}' is not a directory;".format(search_dir))
            sys.exit(1)

    # output files
    samples_fastq_long_tsv = 'samples.fastq.tsv'
    samples_fastq_long_json = 'samples.fastq.json'
    samples_analysis_json = 'samples.analysis.json'
    samples_analysis_tsv = 'samples.analysis.tsv'
    if output_prefix:
        samples_fastq_long_tsv = '{0}.{1}'.format(str(output_prefix), samples_fastq_long_tsv)
        samples_fastq_long_json = '{0}.{1}'.format(str(output_prefix), samples_fastq_long_json)
        samples_analysis_json = '{0}.{1}'.format(str(output_prefix), samples_analysis_json)
        samples_analysis_tsv = '{0}.{1}'.format(str(output_prefix), samples_analysis_tsv)

    # find the R1 fastq files
    fastqs_R1 = []
    for search_dir in search_dirs:
        for fastq_R1 in sorted(find.find(search_dir = search_dir,
                                        inclusion_patterns = [ '*_R1_0*.fastq.gz' ],
                                        search_type = 'file' )):
            fastqs_R1.append(fastq_R1)
    fastqs_R1 = list(sorted(set((fastqs_R1))))

    # parse the files into samples
    samples = []
    for R1_name in fastqs_R1:
        # generate R2 filename
        R2_name = os.path.join(os.path.dirname(R1_name), re.sub(r'(.*)_R1_0([0-9]+\.fastq\.gz)', r'\1_R2_0\2', os.path.basename(R1_name)))
        if not os.path.exists(R2_name): R2_name = None

        # extract sample name
        sample_name = re.sub(r'_S[0-9]{1,3}_L00[0-9]_R1.*', '', os.path.basename(R1_name))
        sample_dict = {
            str(sample_colname): sample_name,
            str(r1_colname): R1_name,
            str(r2_colname): R2_name
        }
        samples.append(sample_dict)

    # save long version of the table; one line per R1 R2 pair
    # with open(samples_fastq_long_tsv, 'w') as f:
    #     writer = csv.DictWriter(f, delimiter= '\t', fieldnames=[str(sample_colname), str(r1_colname), str(r2_colname)])
    #     writer.writeheader()
    #     for item in samples:
    #         writer.writerow(item)

    # save a JSON
    # with open(samples_fastq_long_json, 'w') as f:
    #     json.dump(samples, f, sort_keys = True, indent = 4)

    # reduce to condensed version; one entry per sample with all R1 and R2
    samples_collapsed = collapse.collapse(dicts = samples, collapse_key = str(sample_colname))

    # add some extra metadata
    for sample_dict in samples_collapsed:
        sample_dict[str(tumor_colname)] = sample_dict[str(sample_colname)]
        sample_dict[str(normal_colname)] = str(NA_value)

    # save a JSON
    # with open(samples_analysis_json, 'w') as f:
    #     json.dump(samples_collapsed, f, sort_keys = True, indent = 4)

    # prepare dicts for .tsv printing
    samples_to_print = []
    for sample_dict in samples_collapsed:
        # convert fastq lists to comma separated
        new_R1 = ','.join([str(x) for x in sample_dict[str(r1_colname)]])
        new_R2 = ','.join([str(x) for x in sample_dict['R2']])
        sample_dict[str(r1_colname)] = new_R1
        sample_dict[str(r2_colname)] = new_R2
        samples_to_print.append(sample_dict)

    # write to file
    with open(samples_analysis_tsv, 'w') as f:
        writer = csv.DictWriter(f, delimiter= '\t', fieldnames=[str(sample_colname), str(tumor_colname), str(normal_colname), str(r1_colname), str(r2_colname)])
        writer.writeheader()
        for item in samples_to_print:
            writer.writerow(item)

def parse():
    """
    Parses script arguments
    """
    parser = argparse.ArgumentParser(description='This script will generate samplesheets for the analysis based on .fastq.gz files in the supplied directories')
    parser.add_argument("search_dirs", help="Paths to input samplesheet file", nargs="+")
    parser.add_argument("-p", default = None, dest = 'output_prefix', metavar = 'prefix', help="Prefix for samplesheet files")

    args = parser.parse_args()
    search_dirs = args.search_dirs
    output_prefix = args.output_prefix
    main(search_dirs = search_dirs, output_prefix = output_prefix)

if __name__ == '__main__':
    parse()
