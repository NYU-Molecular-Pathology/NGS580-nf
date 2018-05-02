#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import sqlite3
import csv
import os

# ~~~~~ FUNCTIONS ~~~~~ #
def create_table(conn, table_name, col_name, col_type = "TEXT", is_primary_key = False):
    """
    Create a table in the SQLite if it doesnt exist, with a starting column (required)

    Parameters
    ----------
    conn: sqlite3.Connection object
        connection object to the database
    table_name: str
        the name of the table to create
    col_name: str
        the name of the first column to create in the table
    col_type: str
        the SQLite data type for the column
    is_primary_key: bool
        whether or not the column is the primary key for the table
    """
    with conn:
        cursor = conn.cursor()
        sql_cmd = 'CREATE TABLE IF NOT EXISTS {0}'.format(table_name)
        if is_primary_key:
            table_cmd = ' ({0} {1} PRIMARY KEY)'.format(col_name, col_type)
        else:
            table_cmd = ' ({0} {1})'.format(col_name, col_type)
        sql_cmd = sql_cmd + table_cmd
        print(sql_cmd)
        cursor.execute(sql_cmd)

def sanitize_dict_keys(d):
    """
    Cleans a dictionary's keys for use in the SQLite database as a header by creating a new dict with 'cleaned' keys
    """
    new_dict = { sanitize_str(key): value for key, value in d.items() }
    return(new_dict)

def sanitize_str(string):
    """
    Cleans a character string for use in the database as a header
    """
    string = string.strip().replace(' ', '_')
    string = string.replace('.', '_')
    string = string.replace(':', '_')
    string = string.replace('#', '')
    string = string.replace('%', 'pcnt')
    string = string.replace('/', '_')

    return(string)

def add_column(conn, table_name, col_name, col_type = "TEXT", default_val = None):
    """
    Adds a column to a table

    Parameters
    ----------
    conn: sqlite3.Connection object
        connection object to the database
    table_name: str
        the name of the table in which to create the column
    col_name: str
        the name of the column to create
    col_type: str
        the SQLite data type for the column
    default_val: str
        a default value to use for the column
    """
    sql_cmd = "ALTER TABLE {0} ADD COLUMN '{1}' {2}".format(table_name, col_name, col_type)
    if default_val:
        default_val_cmd = " DEFAULT '{0}'".format(default_val)
        sql_cmd = sql_cmd + default_val_cmd
    try:
        with conn:
            cursor = conn.cursor()
            cursor.execute(sql_cmd)
    except:
        # the column already exists...
        pass

def sqlite_insert(conn, table_name, row, ignore = False, update = False, add_missing_cols = False):
    """
    Inserts a row into a table

    Parameters
    ----------
    conn: sqlite3.Connection object
        connection object to the database
    table_name: str
        the name of the table in which to insert the row
    row: dict
        a dictionary of key: value pairs corresponding to the column names and values of the items in the row to be added
    ignore: bool
        whether the entry should be ignored if it already exists in the table
    add_missing_cols: bool
        whether missing columns should be added to the table. Note: default column type will be used.

    Examples
    --------
    Example usage::

        row = {'key': key, 'value': val}
        sqlite_insert(conn = conn, table_name = vals_table_name, row = row, ignore = True)

    """
    if add_missing_cols:
        colnames = get_colnames(conn = conn, table_name = table_name)
        for key in row.keys():
            if key not in colnames:
                add_column(conn = conn, table_name = table_name, col_name = key)
    cols = ', '.join('"{0}"'.format(col) for col in row.keys())
    vals = ', '.join(':{0}'.format(col) for col in row.keys())
    sql = 'INSERT '
    if ignore:
        sql = sql + 'OR IGNORE '
    sql = sql + 'INTO "{0}" ({1}) VALUES ({2})'.format(table_name, cols, vals)
    with conn:
        conn.cursor().execute(sql, row)

def init_trace_table(trace_file, conn):
    """
    Initializes the 'trace' SQLite table
    """
    with open(trace_file) as f:
        reader = csv.DictReader(f, delimiter = '\t')
        # create the first table
        create_table(conn = conn, table_name = 'trace', col_name = 'hash', col_type = "TEXT", is_primary_key = True)
        for row in reader:
            # re-build dict with clean colname keys
            row = sanitize_dict_keys(d = row)
            # add missing columns to db table
            for key in row.keys():
                add_column(conn = conn, table_name = 'trace', col_name = key, col_type = "TEXT")
            # add the entry to the db
            sqlite_insert(conn = conn, table_name = 'trace', row = row)

def get_hashes(trace_file):
    """
    Gets the hashes from the file
    """
    hashes = []
    with open(trace_file) as f:
        reader = csv.DictReader(f, delimiter = '\t')
        for row in reader:
            hashes.append(row['hash'])
    return(hashes)

def find_hashdir(search_dir, hash):
    """
    """
    hash_dir, hash_subdir = hash.split('/')
    for root, dirs, files in os.walk(search_dir):
        for name in dirs:
            if name == hash_dir:
                hash_dir_path = os.path.join(root, name)
                for sub_root, sub_dirs, sub_files in os.walk(hash_dir_path):
                    for sub_dir in sub_dirs:
                        if sub_dir.startswith(hash_subdir):
                            hash_subdir_path = os.path.join(hash_dir_path, sub_dir)
                            if os.path.exists(hash_subdir_path):
                                return(hash_subdir_path)

def find_nextflow_files(search_dir, nextflow_filenames):
    """
    """
    file_paths = []
    for root, dirs, files in os.walk(search_dir):
        for name in files:
            if name in nextflow_filenames:
                file_path = os.path.join(search_dir, name)
                file_paths.append(file_path)
    return(file_paths)

def get_blob(input_file):
    """
    """
    with open(input_file, 'rb') as f:
        ablob = f.read()
    return(ablob)


# ~~~~~ RUN ~~~~~ #
nextflow_filenames = (
".command.begin",
".command.err",
".command.log",
".command.out",
".command.run",
".command.sh",
".command.stub",
".command.trace",
".exitcode"
)
trace_file = "trace.txt"
sqlite_file = "trace.sqlite"
work_dir = "work"
conn = sqlite3.connect(sqlite_file)
init_trace_table(trace_file, conn)

hashes = { h: { 'path': find_hashdir(search_dir = work_dir, hash = h) } for h in get_hashes(trace_file) }

for hash_key, hash_items in hashes.items():
    hash_items['files'] = find_nextflow_files(search_dir = hash_items['path'], nextflow_filenames = nextflow_filenames)

for hash_key, hash_items in hashes.items():
    hash_items['entry'] = sanitize_dict_keys({ f: get_blob(input_file = f) for f in hash_items['files']})
    hash_items['entry']['hash'] = hash_key
    print(hash_key, hash_items)
    create_table(conn = conn, table_name = sanitize_str(hash_key), col_name = 'hash', col_type = "TEXT", is_primary_key = True)
    for key in hash_items['entry'].keys():
        if key != 'hash':
            add_column(conn = conn, table_name = sanitize_str(hash_key), col_name = key, col_type = "BLOB")
    sqlite_insert(conn = conn, table_name = sanitize_str(hash_key), row = hash_items['entry'])
