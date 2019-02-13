#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
module for .bed file handling
"""
import sys

class Bed(object):
    """

    from bed import Bed; x = Bed("targets.bed"); x.breadthOfCoverage()
    from bed import Bed; x = Bed("targets.bed"); x.hasStrands()
    """
    def __init__(self, path):
        self.path = path

    def hasStrands(self):
        """
        Check if the .bed file has strand information

        $ head targets.bed
        chr1    2985823 2985860 472_145888_63976(PRDM16)_1      0       -
        chr1    3102688 3103038 472_145889_63976(PRDM16)_2      0       -
        """
        with open(self.path) as f:
            for line in f:
                parts = line.split('\t')
                # need at least 6 columns
                if len(parts) < 6:
                    return(False)
                # if has + or - in 6th column, it has strands
                if parts[5].strip() == "+" or parts[5].strip() == "-":
                    return(True)
        return(False)


    def breadthOfCoverage(self):
        """
        Count the total bases covered in the .bed file
        """
        cov = 0
        with open(self.path) as f:
            for line in f:
                chrom, start, stop = line.split('\t')[0:3]
                start = int(start)
                stop = int(stop)
                cov += stop - start
        return(cov)



if __name__ == '__main__':
    args = sys.argv[1:]
    bedFile = args[0]
    command = args[1]
    bedObj = Bed(path = bedFile)
    if command == "hasStrands":
        print(bedObj.hasStrands())
    if command == "breadthOfCoverage":
        print(bedObj.breadthOfCoverage())
