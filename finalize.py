#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script will 'finalize' a Nextflow 'work' directory
by parsing the 'trace.txt' file
identifying subdirs included in the current pipeline
and replacing all files with file stubs
and deleting all other contents
"""
import os
import csv
from bin import find
from multiprocessing import Pool
import argparse

# ~~~~~ CONFIGS ~~~~~ #
configs = {}
configs['work_dir'] = "work"

# ~~~~~~ FUNCTIONS ~~~~~ #
def del_touch_file(fname):
    """
    Destroys the contents a file, leaving an empty file behind at the same path

    Parameters
    ----------
    fname: str
        the path to a file to be destroyed
    """
    if os.path.isfile(fname):
        open(fname, 'wb').close()
    else:
        print("ERROR: item does not exist or is not a file: {0}".format(fname))

def del_touch_worksubdir(subdir_path, exclusion_patterns = ('.command*', '.exitcode')):
    """
    Destroys the contents of all matching files in the subdir. Replaces all files with empty versions of the original

    Parameters
    ----------
    subdir_path: str
        path to the Nextflow 'work' subdirectory to operate on
    exclusion_patterns: list
        iterable of search file patterns to exclude from processing, defaults to ``('.command*', '.exitcode')``
    """
    for item in find.find_gen(search_dir = subdir_path, exclusion_patterns = exclusion_patterns, search_type = 'file', level_limit = 1):
        if not os.path.islink(item):
            print(item)
            del_touch_file(item)

def get_trace_entries(trace_file):
    """
    Loads the trace file entries

    Parameters
    ----------
    trace_file: str
        path to the trace file

    Returns
    -------
    list
        a list of dictionary entries containing the entries from the file
    """
    entries = []
    # load the entries from the trace file
    with open(trace_file) as f:
        reader = csv.DictReader(f, delimiter = '\t')
        for row in reader:
            entries.append(row)
    return(entries)

def get_trace_subdirs(subdir, subsubdir_prefix, work_dir = configs['work_dir']):
    """
    Gets the subdirectories that match the given patterns in the work dir

    Parameters
    ----------
    subdir: str
        the path to the subdir in the work directory, e.g. ``work/b8``
    subsubdir_prefix: str
        the hash pattern of the sub-subdir to search for, e.g. ``8f0cbc1c``
    work_dir: str
        the Nextflow 'work' directory to process
    """
    subdir_path = os.path.join(work_dir, subdir)  #work/b8
    if not os.path.isdir(subdir_path):
        print("ERROR: path is not a dir or does not exist: {0}".format(subdir_path))
        raise
    subsubdir_pattern = '{0}*'.format(subsubdir_prefix) # 8f0cbc1c*
    # get the sub-subdirs that match
    for item in find.find_gen(search_dir = subdir_path, inclusion_patterns = [subsubdir_pattern], search_type = 'dir', level_limit = 1):
        yield(item) # work/b8/8f0cbc1c7f65e4b5eab2e7be4105ac

def process_subdirs(subdir, subsubdir_prefix):
    """
    Applies the processing functions to all subdirs found that match the given subdir criteria

    Parameters
    ----------
    subdir: str
        the path to the subdir in the work directory, e.g. ``work/b8``
    subsubdir_prefix: str
        the hash pattern of the sub-subdir to search for, e.g. ``8f0cbc1c``

    Returns
    -------
    list
        a list of all matching subdir paths that were processed
    """
    sub_subdirs = []
    for sub_subdir in get_trace_subdirs(subdir = subdir, subsubdir_prefix = subsubdir_prefix):
        del_touch_worksubdir(subdir_path = sub_subdir)
        sub_subdirs.append(sub_subdir)
    return(sub_subdirs)


def main(trace_file):
    """
    Main control function for the script
    """
    print("Parsing trace file: {0}".format(trace_file))
    entries = get_trace_entries(trace_file = trace_file)

    subdirs = []
    pool = Pool(processes=4)


    # parse the entries for dirs to process
    for entry in entries:
        task_hash = entry["hash"] # b8/8f0cbc1c
        parts = task_hash.split('/')
        subdir = parts[0] # b8
        subsubdir_prefix = parts[1] # 8f0cbc1c
        subdirs.append((subdir, subsubdir_prefix))

    # destroy the subdir contents
    print("{0} subdir patterns found".format(len(subdirs)))
    print("Removing file contents in subdirs...")
    results = [pool.apply_async(process_subdirs, args=(subdir, subsubdir_prefix)) for subdir, subsubdir_prefix in subdirs]
    output = [p.get() for p in results]

    # list of all the subdirs that were cleaned up
    processed_subdirs = set([item for sublist in output for item in sublist])
    print("{0} subdirs were processed".format(len(processed_subdirs)))

    # list of all subdirs not included in the trace file patterns
    unproccessed_subdirs = []
    for item in find.find_gen(search_dir = configs['work_dir'], search_type = 'dir', level_limit = 2):
        if os.path.isdir(item) and not os.path.islink(item) and item not in processed_subdirs and item.count(os.path.sep) == 2:
            unproccessed_subdirs.append(item)
    print("{0} subdirs in the 'work' dir were not processed: ".format(len(unproccessed_subdirs)))
    print(unproccessed_subdirs)

def parse():
    """
    Parses script arguments to run the main function
    """
    parser = argparse.ArgumentParser(description='Destroys the contents of all pipeline files in the Nextflow "work" directory and replaces them with an empty file stub')
    parser.add_argument("-t", "--trace-file", default = "trace-NGS580.txt", dest = 'trace_file', metavar = 'trace_file', help="Nextflow trace file")
    parser.add_argument("-w", "--work-dir", default = "work", dest = 'work_dir', metavar = 'work_dir', help="Nextflow work directory")

    args = parser.parse_args()

    configs['work_dir'] = args.work_dir
    trace_file = args.trace_file

    main(trace_file = trace_file)


# ~~~~~ RUN ~~~~~ #
if __name__ == '__main__':
    parse()
