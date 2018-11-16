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

    - ``samples.fastq.tsv``: one line per R1 R2 fastq file pair (optional)
    - ``samples.fastq.json``: one entry per R1 R2 fastq file pair (optional)
    - ``samples.analysis.tsv``: one line per sample with R1 R2 file pairs and metadata (default)
    - ``samples.analysis.json``: one entry per sample with R1 R2 file pairs and metadata (optional)

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
def get_samplename(fastq, mode = "noLaneSplit"):
    """
    Derives a sample name from a .fastq filename, using Illumina standard filenaming syntax

    Parameters
    ----------
    fastq: str
        a character string representing the path to a file to be parsed
    mode: str
        one of 'LaneSplit' or 'noLaneSplit'; whether or not the fastq files were produced with lane splitted enabled in bcl2fastq. A fastq with lane splitting will have a filename such as 'SampleName_S1_L001_R1_001.fastq.gz' while a fastq produced without lane splitting will be named such as 'SampleName_S1_R1_001.fastq.gz'

    Returns
    -------
    str
        A string representing the derived sample name
    """
    if mode == "LaneSplit":
        sample_name = re.sub(r'_S[0-9]{1,3}_L00[0-9]_R[1-2].*', '', os.path.basename(fastq))
    elif mode == "noLaneSplit":
        sample_name = re.sub(r'_S[0-9]{1,3}_R[1-2].*', '', os.path.basename(fastq))
    else:
        print("ERROR: invalid mode")
        raise
    return(sample_name)

def write_tsv(samples, output_file, fieldnames, append = False):
    """
    save a TSV file

    Parameters
    ----------
    samples: list
        a list of dictionaries to be written to the file
    output_file: str
        path to output file
    fieldnames: list
        list of column headers to use for the output file1
    append: bool
        ``True`` or ``False``, whether to append to the current output file
    """
    if append:
        fout = open(output_file, "a")
        writer = csv.DictWriter(fout, delimiter= '\t', fieldnames = fieldnames)
    else:
        fout = open(output_file, "w")
        writer = csv.DictWriter(fout, delimiter= '\t', fieldnames = fieldnames)
        writer.writeheader()
    for item in samples:
        writer.writerow(item)
    fout.close()

def write_json(samples, output_file):
    """
    save a JSON file

    Parameters
    ----------
    samples: list
        a list of dictionaries to be written to the file
    output_file: str
        path to output file
    """
    with open(output_file, 'w') as f:
        json.dump(samples, f, sort_keys = True, indent = 4)

def main(**kwargs):
    """
    Main control function for the script
    """
    # get args
    search_dirs = kwargs.pop('search_dirs')
    output_prefix = kwargs.pop('output_prefix', None)
    NA_value = kwargs.pop('NA_value', "NA")
    tumor_colname = kwargs.pop('tumor_colname', 'Tumor')
    normal_colname = kwargs.pop('normal_colname', 'Normal')
    r1_colname = kwargs.pop('r1_colname', 'R1')
    r2_colname = kwargs.pop('r2_colname', 'R2')
    sample_colname = kwargs.pop('sample_colname', 'Sample')
    name_mode = kwargs.pop('name_mode', 'noLaneSplit')
    write_long_tsv = kwargs.pop('write_long_tsv', False)
    write_long_json = kwargs.pop('write_long_json', False)
    write_analysis_json = kwargs.pop('write_analysis_json', False)
    samples_fastq_long_tsv = kwargs.pop('samples_fastq_long_tsv', 'samples.fastq.tsv')
    samples_fastq_long_json = kwargs.pop('samples_fastq_long_json', 'samples.fastq.json')
    samples_analysis_json = kwargs.pop('samples_analysis_json', 'samples.analysis.json')
    samples_analysis_tsv = kwargs.pop('samples_analysis_tsv', 'samples.analysis.tsv')
    append = kwargs.pop('append', False)

    # validate inputs
    if len(search_dirs) < 1:
        print("ERROR: no directories provided")
        sys.exit(1)
    for search_dir in search_dirs:
        if not os.path.isdir(search_dir):
            print("ERROR: '{0}' is not a directory;".format(search_dir))
            sys.exit(1)

    if output_prefix:
        samples_fastq_long_tsv = '{0}.{1}'.format(str(output_prefix), samples_fastq_long_tsv)
        samples_fastq_long_json = '{0}.{1}'.format(str(output_prefix), samples_fastq_long_json)
        samples_analysis_json = '{0}.{1}'.format(str(output_prefix), samples_analysis_json)
        samples_analysis_tsv = '{0}.{1}'.format(str(output_prefix), samples_analysis_tsv)

    # find the R1 fastq files
    # TODO: clean this up
    fastqs_R1 = []
    for search_dir in search_dirs:
        for fastq_R1 in sorted(find.find(search_dir = search_dir,
                                        inclusion_patterns = [ '*_R1_0*.fastq.gz' ],
                                        search_type = 'file' )):
            fastqs_R1.append(fastq_R1)
    fastqs_R1 = list(sorted(set((fastqs_R1))))

    # parse the files into samples
    # TODO: clean this up too
    samples = []
    for R1_name in fastqs_R1:
        # generate R2 filename
        R2_name = os.path.join(os.path.dirname(R1_name), re.sub(r'(.*)_R1_0([0-9]+\.fastq\.gz)', r'\1_R2_0\2', os.path.basename(R1_name)))
        if not os.path.exists(R2_name): R2_name = None

        # extract sample name
        sample_name = get_samplename(fastq = R1_name, mode = name_mode)
        sample_dict = {
            str(sample_colname): sample_name,
            str(r1_colname): R1_name,
            str(r2_colname): R2_name
        }
        samples.append(sample_dict)

    # save long version of the table; one line per R1 R2 pair
    if write_long_tsv:
        write_tsv(samples = samples, output_file = samples_fastq_long_tsv, fieldnames=[str(sample_colname), str(r1_colname), str(r2_colname)], append = append)

    # save a JSON
    if write_long_json:
        write_json(samples = samples, output_file = samples_fastq_long_json)

    # reduce to condensed version; one entry per sample with all R1 and R2
    samples_collapsed = collapse.collapse(dicts = samples, collapse_key = str(sample_colname))

    # add some extra metadata
    for sample_dict in samples_collapsed:
        sample_dict[str(tumor_colname)] = sample_dict[str(sample_colname)]
        sample_dict[str(normal_colname)] = str(NA_value)

    # save a JSON
    if write_analysis_json:
        write_json(samples = samples_collapsed, output_file = samples_analysis_json)

    # prepare dicts for .tsv printing
    samples_to_print = []
    for sample_dict in samples_collapsed:
        # convert fastq lists to comma separated
        new_R1 = ','.join([str(x) for x in sample_dict[str(r1_colname)]])
        new_R2 = ','.join([str(x) for x in sample_dict['R2']])
        sample_dict[str(r1_colname)] = new_R1
        sample_dict[str(r2_colname)] = new_R2
        samples_to_print.append(sample_dict)

    # write to file; `samples.analysis.tsv`
    write_tsv(samples = samples_to_print, output_file = samples_analysis_tsv, fieldnames=[str(sample_colname), str(tumor_colname), str(normal_colname), str(r1_colname), str(r2_colname)], append = append)

def parse():
    """
    Parses script arguments
    """
    parser = argparse.ArgumentParser(description='This script will generate samplesheets for the analysis based on .fastq.gz files in the supplied directories')
    parser.add_argument("search_dirs", help="Paths to input directories to search for .fastq files to use for samplesheet creation", nargs="+")
    parser.add_argument("-p", default = None, dest = 'output_prefix', metavar = 'prefix', help="Prefix for samplesheet files")
    parser.add_argument("--name-mode", default = 'noLaneSplit', dest = 'name_mode', metavar = 'filename mode', help="Mode for parsing fastq filenames. Default: 'noLaneSplit', alternative: 'LaneSplit'")
    parser.add_argument("--append", action = 'store_true', dest = 'append', help="Append newly discovered samples to existing samplesheet")

    args = parser.parse_args()
    main(**vars(args))

if __name__ == '__main__':
    parse()
