#!/usr/bin/env python
# -*- coding: utf-8 -*-
def collapse(dicts, collapse_key):
    """
    Collapses a list of dicts based on a given key

    Parameters
    ----------
    dicts: list
        a list of dicts that are assumed to have common keys
    collapse_key: str
        dict key to collapse on

    Returns
    -------
    list
        a list of dicts

    Examples
    --------
    Example usage

        x = [{'ID': 'foo', 'File1': 'foo-1-1.txt', 'File2': 'foo-1-2.txt'},
        {'ID': 'foo', 'File1': 'foo-2-1.txt', 'File2': 'foo-2-2.txt'},
        {'ID': 'foo', 'File1': 'foo-3-1.txt', 'File2': 'foo-3-2.txt'},
        {'ID': 'bar', 'File1': 'bar-1-1.txt', 'File2': 'bar-1-2.txt'},
        {'ID': 'bar', 'File1': 'bar-2-1.txt', 'File2': 'bar-2-2.txt'},
        {'ID': 'bar', 'File1': 'bar-3-1.txt', 'File2': 'bar-3-2.txt'}]
        collapse(x, 'ID')

        # {'File2': ['foo-1-2.txt', 'foo-2-2.txt', 'foo-3-2.txt'], 'File1': ['foo-1-1.txt', 'foo-2-1.txt', 'foo-3-1.txt'], 'ID': 'foo'}
        # {'File2': ['bar-1-2.txt', 'bar-2-2.txt', 'bar-3-2.txt'], 'File1': ['bar-1-1.txt', 'bar-2-1.txt', 'bar-3-1.txt'], 'ID': 'bar'}
    """
    ids = list(set([i[collapse_key] for i in dicts])) # unique values from all dicts for the given key
    # keys = set([k for k in i.keys() for i in dicts]) # all keys in all dicts
    # for some reason this stopped working... ? something to do with Python 3.6??
    keys = []
    for d in dicts:
        for k in d.keys():
            keys.append(k)
    keys = list(set(keys))
    collapsed_list = [] # list to hold output dicts
    for i in ids:
        y = {collapse_key: i} # retain value of collapsed field
        for k in keys :
            if k != collapse_key: # add the rest of the keys
                if k not in y.keys():
                    y.update({k: []}) # initialize them as lists for appending
                for z in dicts: # iterate over the original list of dicts again
                    if z[collapse_key] == i: # only the entries that match the selected collapse value
                        y[k].append(z[k]) # append value to the list
        collapsed_list.append(y) # add the collapsed dict to the list
    return(collapsed_list)
