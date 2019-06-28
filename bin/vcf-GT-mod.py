#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Add GT field to .vcf file, in order to fix issues with programs that require this field

e.g.: http://seqanswers.com/forums/showthread.php?t=44799

Expects a .vcf file with 'TUMOR' and 'NORMAL' sample labels
"""
import sys

def main(input_vcf, output_vcf):
    """
    main control function for the script
    """
    # header line to insert
    genotype_FORMAT_str = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'

    # read input vcf lines
    input_lines = []
    with open(input_vcf) as fin:
        for line in fin:
            input_lines.append(line)

    if len(input_lines) < 1:
        print("ERROR: no input lines read from file: {0}".format(input_vcf))
        raise

    # find index of first FORMAT line in header
    first_FORMAT_index = 0
    for i, line in enumerate(input_lines):
        if line.startswith('##FORMAT='):
            break
        else:
            first_FORMAT_index += 1

    # insert the new header line
    input_lines.insert(first_FORMAT_index, genotype_FORMAT_str)

    # find the vcf field header line;
    # '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOR\n'
    vcf_fields_index = 0
    for i, line in enumerate(input_lines):
        if line.startswith('#CHROM'):
            break
        else:
            vcf_fields_index += 1

    # get the vcf fields headers line
    vcf_fields_header = input_lines[vcf_fields_index]

    # figure out positions of TUMOR and NORMAL columns
    vcf_fields = vcf_fields_header.split('\t')
    tumor_index = [ n.strip() for n in vcf_fields ].index('TUMOR')
    normal_index = [ n.strip() for n in vcf_fields ].index('NORMAL')
    format_index = [ n.strip() for n in vcf_fields ].index('FORMAT')

    # update the lines for the file;
    # old:
    # 'chr21	42871545	.	A	AT	.	LowEVS	SOMATIC;QSI=69;TQSI=1;NT=het;QSI_NT=69;TQSI_NT=1;SGT=het->ref;MQ=60;MQ0=0;RU=T;RC=5;IC=6;IHP=5;SomaticEVS=0	DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50:BCN50	48:48:21,22:27,27:0,0:48.44:0.65:0:0	56:56:47,47:0,0:9,9:61.98:1.08:0:0.02\n'
    # new:
    # 'chr21	42871545	.	A	AT	.	LowEVS	SOMATIC;QSI=69;TQSI=1;NT=het;QSI_NT=69;TQSI_NT=1;SGT=het->ref;MQ=60;MQ0=0;RU=T;RC=5;IC=6;IHP=5;SomaticEVS=0	GT:DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50:BCN50	0/0:48:48:21,22:27,27:0,0:48.44:0.65:0:0	0/1:56:56:47,47:0,0:9,9:61.98:1.08:0:0.02\n'
    new_lines = []
    for line in input_lines:
        if line.startswith('#'):
            # skip header lines
            new_lines.append(line)
        else:
            # modify variant lines
            old_parts = line.split('\t')
            new_parts = []
            for i, part in enumerate(old_parts):
                new_part = part
                if i == format_index:
                    # prepend 'GT:' to the FORMAT field
                    new_part = 'GT:' + part
                elif i == normal_index:
                    # prepend '0/0:' for NORMAL
                    new_part = '0/0:' + part
                elif i == tumor_index:
                    # prepend '0/1:' for TUMOR
                    new_part = '0/1:' + part
                new_parts.append(new_part)
            new_line = '\t'.join(new_parts)
            new_lines.append(new_line)

    # write to new file
    with open(output_vcf, "w") as fout:
        for line in new_lines:
            fout.write(line)

if __name__ == '__main__':
    args = sys.argv[1:]
    input_vcf = args[0]
    output_vcf = args[1]
    main(input_vcf, output_vcf)
