#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Retrieves entries from a .bed file that match the provided chromosome.
Prints each matching line to stdout

Example usage:
$ ./subset_bed.py chr1 ../regions.bed
"""
import sys

def main(chrom, bed_file):
    """
    Main control function for the script
    """
    with open(bed_file) as f:
        for line in f:
            if len(line.split()) > 0:
                if line.split()[0] == chrom:
                    sys.stdout.write(line)
def parse():
    """
    Parse script args if called as script
    """
    chrom = sys.argv[1]
    bed_file = sys.argv[2]
    main(chrom, bed_file)

if __name__ == '__main__':
    parse()
