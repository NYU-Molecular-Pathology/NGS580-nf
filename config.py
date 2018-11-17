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

    for key, new_value in data.items():
        if key in old_data:
            # only update non-None values
            if new_value is not None:
                # if the old data value was a list, append new items
                if isinstance(old_data[key], (list,)):
                    # start a new list from the old items
                    new_list = [ item for item in old_data[key] ]
                    # check if new item is single value or list of values
                    if isinstance(new_value, (list,)):
                        # add each new value
                        for item in new_value:
                            new_list.append(item)
                    else:
                        # add the new value
                        new_list.append(new_value)
                    old_data[key] = list(set(new_list))
                # otherwise, overwrite old value
                else:
                    old_data[key] = new_value

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
    updateFile = kwargs.pop('updateFile', False)

    if fastqDirs is None:
        fastqDirs = []

    data = {
    'runID': runID,
    'fastqDirs': list(set(fastqDirs)),
    'samplesheet': samplesheet
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
    parser.add_argument("-u", "--update", default = False, dest = 'updateFile', help="JSON file to update")

    args = parser.parse_args()

    main(**vars(args))

if __name__ == '__main__':
    parse()
