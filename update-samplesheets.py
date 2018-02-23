#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script will update the samplesheets ``samples.analysis.json`` and ``samples.analysis.tsv`` output by the ``generate-samplesheets.py`` script.
Use this script to update those files with tumor-normal pairs sample metadata, from the file ``samples.tumor.normal.csv`` (output by Excel)


"""
import csv
import json

# ~~~~~ CONFIGS ~~~~~ #
NA_value = "NA"
tumor_normal_samples_delim = ','
tumor_normal_sheet = 'samples.tumor.normal.csv'
samples_analysis_json = 'samples.analysis.json'
samples_analysis_tsv = 'samples.analysis.tsv'
tsv_delim = '\t'
NA_value = "NA"
tumor_colname = 'Tumor'
normal_colname = 'Normal'
r1_colname = 'R1'
r2_colname = 'R2'
sample_colname = 'Sample'

# ~~~~~ FUNCTIONS ~~~~~ #
def update_samples(old_data, tumor_normal_samples):
    """
    Updates the old data read in from the analysis files with the new tumor normal samples values

    Parameters
    ----------
    old_data: list
        a list of dicts with the values read in from the analysis files
    tumor_normal_samples: list
        a list of dicts with the values read in from the tumor normal pairs sheet

    Returns
    -------
    list
        a list of dicts with the updated data
    """
    # parse the input samples
    for item in tumor_normal_samples:
        tumor_ID = item[tumor_colname]
        normal_ID = item[normal_colname]

        for dat in old_data:
            # find samples that match the new tumor ID
            if dat[sample_colname] == tumor_ID:
                # update the tumor and normal values
                dat[tumor_colname] = tumor_ID
                dat[normal_colname] = normal_ID
    return(old_data)

def update_analysis_json(input_json, tumor_normal_samples, overwrite = True, output_json = 'new.json'):
    """
    Updates the 'Normal' value in the ``input_json`` for all matching entries from ``tumor_normal_samples``

    Parameters
    ----------
    input_json: str
        path to input JSON file to be updated
    output_json: str
        path to save JSON output
    tumor_normal_samples: list
        list of dicts with the values to update in the JSON, must match the input JSON format
    overwrite: bool
        whether to overwrite the old file with the new one
    """
    if overwrite:
        output_json = input_json

    # load data from JSON
    with open(input_json) as data_file:
        data = json.load(data_file)

    data = update_samples(old_data = data, tumor_normal_samples = tumor_normal_samples)

    # save the output
    with open(output_json, 'w') as f:
        json.dump(data, f, sort_keys = True, indent = 4)


def update_analysis_tsv(input_tsv, tumor_normal_samples, overwrite = True, output_tsv = 'new.tsv', input_delim = '\t'):
    """
    """
    if overwrite:
        output_tsv = input_tsv

    # load data from .TSV
    data = []
    with open(input_tsv) as f:
        reader = csv.DictReader(f, delimiter = input_delim)
        for row in reader:
            data.append(row)

    data = update_samples(old_data = data, tumor_normal_samples = tumor_normal_samples)

    # save the output
    output_fields = [sample_colname, tumor_colname, normal_colname, r1_colname, r2_colname]
    with open(output_tsv, 'w') as f:
        writer = csv.DictWriter(f, delimiter = '\t', fieldnames = output_fields)
        writer.writeheader()
        for item in data:
            writer.writerow(item)

def main():
    """
    Main control function for the script
    """
    # load input sheet to update with
    tumor_normal_samples = []
    with open(tumor_normal_sheet) as f:
        reader = csv.DictReader(f, delimiter = tumor_normal_samples_delim)
        for row in reader:
            sample_dict = {tumor_colname: row[tumor_colname], normal_colname: row[normal_colname]}
            tumor_normal_samples.append(sample_dict)

    # update_analysis_json(input_json = samples_analysis_json, tumor_normal_samples = tumor_normal_samples)
    update_analysis_tsv(input_tsv = samples_analysis_tsv, tumor_normal_samples = tumor_normal_samples)

if __name__ == '__main__':
    main()
