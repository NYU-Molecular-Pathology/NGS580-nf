#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Creates an updated config file for the pipeline
"""
import json
import argparse

def make_config(data):
    """
    Prints JSON config to stdout
    """
    print(json.dumps(data, indent=4))

def update_config(updateFile, data):
    """
    Update the values of a JSON config file for non-null values in the provided data
    """
    f = open(updateFile)
    old_data = json.load(f)
    f.close()

    for key, value in data.items():
        if key in old_data:
            if value is not None:
                old_data[key] = value

    f = open(updateFile, "w")
    json.dump(old_data, f, indent = 4)
    f.close()

def main(**kwargs):
    """
    Main control function for the script
    """
    runID = kwargs.pop('runID', None)
    fastqDirs = kwargs.pop('fastqDirs', [])
    samplesheet = kwargs.pop('samplesheet', None)
    pairsheet = kwargs.pop('pairsheet', None)
    updateFile = kwargs.pop('updateFile', False)

    if fastqDirs is None:
        fastqDirs = []

    data = {
    'runID': runID,
    'fastqDirs': fastqDirs,
    'samplesheet': samplesheet,
    'pairsheet': pairsheet
    }

    if updateFile is False:
        make_config(data)
    else:
        update_config(updateFile, data)

def parse():
    """
    Parses script args
    """
    parser = argparse.ArgumentParser(description='Creates an updated config file for the pipeline')
    parser.add_argument("--runID", default = None, dest = 'runID', help="Run ID")
    parser.add_argument("--fastqDirs", nargs = '+', dest = 'fastqDirs', help="Run directory")
    parser.add_argument("--samplesheet", default = None, dest = 'samplesheet', help="Samplesheet file")
    parser.add_argument("--pairsheet", default = None, dest = 'pairsheet', help="Tumor/normal pairs samplesheet file")
    parser.add_argument("-u", "--update", default = False, dest = 'updateFile', help="JSON file to update")

    args = parser.parse_args()

    main(**vars(args))

if __name__ == '__main__':
    parse()
