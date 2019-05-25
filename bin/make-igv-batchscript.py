#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script to create an IGV batchscript
https://software.broadinstitute.org/software/igv/PortCommands
http://software.broadinstitute.org/software/igv/automation
https://software.broadinstitute.org/software/igv/batch

example IGV batch script:

new
snapshotDirectory IGV_Snapshots
load test_alignments.bam
genome hg19
maxPanelHeight 500
goto chr1:713167-714758
snapshot chr1_713167_714758_h500.png
goto chr1:713500-714900
snapshot chr1_713500_714900_h500.png
exit

Usage:
./make-batchscript.py foo.bam bar.bam

"""
import os
import argparse

def append_string(string, output_file):
    """
    Append a string to a file
    """
    with open(output_file, "a") as myfile:
        myfile.write(string + '\n')

def make_regions(regions_file):
    """
    Parse the .bed format regions file to generate the IGV location and output filenames
    """
    regions = []
    with open(regions_file) as f:
        for line in f:
            if len(line.split()) >= 3:
                chrom, start, stop = line.split()[0:3]
            elif len(line.split()) == 2:
                chrom, start = line.split()
                stop = start
            # make IGV format location
            loc = '{0}:{1}-{2}'.format(chrom, start, stop)
            filename = '{0}_{1}_{2}.png'.format(chrom, start, stop)
            region = {'chrom': chrom, 'start': start, 'stop': stop, 'loc': loc, 'filename': filename}
            regions.append(region)
    return(regions)

def main(**kwargs):
    """
    Main control function for the script
    """
    input_files = kwargs.pop('input_files')
    regions_file = kwargs.pop('regions_file', "regions.bed")
    snapshotDirectory = kwargs.pop('snapshotDirectory', "IGV_snapshots")
    batchscript_file = kwargs.pop('batchscript_file', "IGV_snapshots.bat")
    image_height = int(kwargs.pop('image_height', 500))
    genome = kwargs.pop('genome', "hg19")

    regions = make_regions(regions_file)

    append_string("new", batchscript_file)
    append_string("snapshotDirectory " + snapshotDirectory, batchscript_file)
    append_string("genome " + genome, batchscript_file)
    for input_file in input_files:
        append_string("load " + input_file, batchscript_file)
    append_string("maxPanelHeight " + str(image_height), batchscript_file)
    for region in regions:
        append_string("goto " + region['loc'], batchscript_file)
        append_string("snapshot " + region['filename'], batchscript_file)
    append_string("exit", batchscript_file)

def parse():
    """
    Parses script args
    """
    parser = argparse.ArgumentParser(description='IGV batchscript creator')
    parser.add_argument("input_files",
        nargs='+',
        help="pathes to the files to create snapshots from e.g. .bam, .bigwig, etc.")
    parser.add_argument("-r", "--regions",
        default = "regions.bed",
        dest = 'regions_file',
        metavar = 'regions_file',
        help="Path to .bed formatted regions file")
    parser.add_argument("-b",
        default = "IGV_snapshots.bat",
        dest = 'batchscript_file',
        metavar = 'batchscript_file',
        help="Name of the IGV batchscript file to create")
    parser.add_argument("-d",
        default = "IGV_snapshots.bat",
        dest = 'snapshotDirectory',
        metavar = 'snapshotDirectory',
        help="Name of the IGV snapshot directory to save images to")
    parser.add_argument("--height",
        default = 500,
        dest = 'image_height',
        metavar = 'image_height',
        help="Height in pixels of the images to create")
    parser.add_argument("--genome",
        default = "hg19",
        dest = 'genome',
        metavar = 'genome',
        help="Name of genome to use in IGV")

    args = parser.parse_args()
    main(**vars(args))

if __name__ == '__main__':
    parse()
