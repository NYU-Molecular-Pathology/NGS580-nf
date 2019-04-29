#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Output a JSON with information about the git repo at the desired location (default: pwd)

Tested with git version 2.16.0
"""
import os
import sys
import subprocess as sp
import argparse
import json

class SubprocessCmd(object):
    """
    A command to be run in subprocess

    Examples
    --------
    Example usage::

        run_cmd = SubprocessCmd(command = ['echo', 'foo']).run()
        print(run_cmd.proc_stdout)
        print(run_cmd.proc_stderr)
    """
    def __init__(self, command):
        self.command = command

    def run(self, command = None):
        """
        Run the command, capture the process object

        Parameters
        ----------
        command: list
            a list of string args to be run in the subprocess

        Notes
        -----
        `universal_newlines=True` required for Python 2 3 compatibility with stdout parsing
        """
        if not command:
            command = self.command
        if command:
            self.process = sp.Popen(command, stdout = sp.PIPE, stderr = sp.PIPE, shell = False, universal_newlines = True)
            self.proc_stdout, self.proc_stderr = self.process.communicate()
            self.proc_stdout = self.proc_stdout.strip()
            self.proc_stderr = self.proc_stderr.strip()
        # else:
        #     logger.error('No command supplied')
        # TODO: what to do here?
        return(self)


class DirHop(object):
    """
    A class for executing commands in the context of a different working directory
    adapted from: https://mklammler.wordpress.com/2011/08/14/safe-directory-hopping-with-python/

    Parameters
    ----------
    directory: str
        path to directory to execute commands from

    Examples
    --------
    with DirHop('/some/dir') as d:
        do_something()
    """
    def __init__(self, directory):
        self.old_dir = os.getcwd()
        self.new_dir = directory
    def __enter__(self):
        os.chdir(self.new_dir)
        return(self)
    def __exit__(self, type, value, traceback):
        os.chdir(self.old_dir)
        return(isinstance(value, OSError))

def main(**kwargs):
    """
    Main control function for the script
    """
    repo_dir = kwargs.pop('repo_dir', '.')
    output_json = kwargs.pop('output_json', 'git.json')
    repo_path = os.path.realpath(repo_dir)

    if output_json == "None" or output_json == "stdout":
        fout = sys.stdout
    else:
        fout = open(output_json, "w")

    data = {}

    current_branch_cmd = ['git', 'rev-parse', '--abbrev-ref', 'HEAD']
    current_commit_cmd = ['git', 'rev-parse', '--short', 'HEAD']
    current_tag_cmd = ['git', 'describe', '--tags']
    recent_tag_cmd = ['git', 'describe', '--tags', '--abbrev=0']

    with DirHop(repo_path) as d:
        current_branch = SubprocessCmd(command = current_branch_cmd).run()
        data['currentBranch'] = current_branch.proc_stdout

        current_commit = SubprocessCmd(command = current_commit_cmd).run()
        data['currentCommit'] = current_commit.proc_stdout

        current_tag = SubprocessCmd(command = current_tag_cmd).run()
        data['currentTag'] = current_tag.proc_stdout

        recent_tag = SubprocessCmd(command = recent_tag_cmd).run()
        data['tag'] = recent_tag.proc_stdout

    json.dump(data, fout, indent = 4)

def parse():
    """
    Parses script args
    """
    parser = argparse.ArgumentParser(description='Saves information about git repo at the given location')
    parser.add_argument("-d", "--dir", default = '.', dest = 'repo_dir', help="Path to directory with git repository")
    parser.add_argument("-o", default = "git.json", dest = 'output_json', help="Output JSON file")
    args = parser.parse_args()

    main(**vars(args))

if __name__ == '__main__':
    parse()
