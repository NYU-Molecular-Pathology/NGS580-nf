#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions for finding files and dirs

adapted from https://github.com/NYU-Molecular-Pathology/util/blob/8bd30b89aa82706f0a10d329fc6a99eb107692a9/find.py
"""

import os
import sys
import itertools
import fnmatch
from collections import defaultdict

def find(search_dir, inclusion_patterns = ('*',), exclusion_patterns = (), search_type = 'all', num_limit = None, level_limit = None, match_mode = "any"):
    """
    Function to search for files and directories

    Parameters
    ----------
    search_dir: str
        path to the directory in which to search for files and subdirectories
    inclusion_patterns: list or tuple
        a list or tuple of patterns to match files/dirs against for inclusion in match output
    exclusion_patterns: list or tuple
        a list or tuple of patterns to match files/dirs against for exclusion from match output
    num_limit: int
        the number of matches to return; use `None` for no limit
    level_limit: int
        the number of directory levels to recurse; 0 is parent dir only
    match_mode:
        'any' or 'all'; matches any of the provided inclusion_patterns, or all of them
    search_type:
        'all', 'file', or 'dir'; type of items to find

    Returns
    -------
    list
        a list of matching file or directory paths
    """
    import sys
    import itertools
    if num_limit != None:
        matches = []
        for item in find_gen(search_dir = search_dir, inclusion_patterns = inclusion_patterns, exclusion_patterns = exclusion_patterns, search_type = search_type, level_limit = level_limit, match_mode = match_mode):
            if len(matches) < int(num_limit):
                matches.append(item)
        return(matches)
    else:
        matches = [item for item in find_gen(search_dir = search_dir, inclusion_patterns = inclusion_patterns, exclusion_patterns = exclusion_patterns, search_type = search_type, level_limit = level_limit, match_mode = match_mode)]
        return(matches)

def find_gen(search_dir, inclusion_patterns = ('*',), exclusion_patterns = (), search_type = 'all', level_limit = None, match_mode = "any"):
    """
    Generator function to return file matches. Used internally by `find`

    Parameters
    ----------
    search_dir: str
        path to the directory in which to search for files and subdirectories
    inclusion_patterns: list or tuple
        a list or tuple of patterns to match files/dirs against for inclusion in match output
    exclusion_patterns: list or tuple
        a list or tuple of patterns to match files/dirs against for exclusion from match output
    level_limit: int
        the number of directory levels to recurse; 0 is parent dir only
    match_mode:
        'any' or 'all'; matches any of the provided inclusion_patterns, or all of them
    search_type:
        'all', 'file', or 'dir'; type of items to find
    """
    import os
    import sys
    import fnmatch
    search_dir = search_dir.rstrip(os.path.sep)
    num_sep = search_dir.count(os.path.sep)
    for root, dirs, files in os.walk(search_dir):
        # choose which items to search
        if search_type == 'all':
            items = dirs + files
        elif search_type == 'dir':
            items = dirs
        elif search_type == 'file':
            items = files
        else:
            print("ERROR: Search type '{0}' not valid, exiting script".format(search_type))
            sys.exit()
        # yeild the results
        for item in super_filter(names = items, inclusion_patterns = inclusion_patterns, exclusion_patterns = exclusion_patterns, match_mode = match_mode):
            yield(os.path.join(root, item))
        # check for a level limit
        if level_limit != None:
            num_sep_this = root.count(os.path.sep)
            if num_sep + int(level_limit) <= num_sep_this:
                del dirs[:]


def super_filter(names, inclusion_patterns = ('*',), exclusion_patterns = (), match_mode = "any"):
    """
    Enhanced version of `fnmatch.filter()` that accepts multiple inclusion and exclusion patterns.

    Filter the input names by choosing only those that are matched by
    some pattern in `inclusion_patterns` _and_ not by any in `exclusion_patterns`.

    Adapted from:
    https://codereview.stackexchange.com/questions/74713/filtering-with-multiple-inclusion-and-exclusion-patterns
    """
    included = multi_filter(names, patterns = inclusion_patterns, match_mode = match_mode)
    excluded = multi_filter(names, patterns = exclusion_patterns, match_mode = match_mode)
    for item in set(included) - set(excluded):
        yield(item)

def multi_filter(names, patterns, match_mode = "any"):
    """
    Generator function which yields the names that match one or more of the patterns.
    """
    for name in names:
        basename = os.path.basename(name)
        # in case a single string was passed as a pattern
        if isinstance(patterns, str):
            if fnmatch.fnmatch(basename, patterns):
                yield(name)
        # patterns is not an empty list
        elif patterns:
            if match_mode == 'any':
                if any(fnmatch.fnmatch(basename, pattern) for pattern in patterns):
                    yield(name)
            elif match_mode == 'all':
                if all(fnmatch.fnmatch(basename, pattern) for pattern in patterns):
                    yield(name)
