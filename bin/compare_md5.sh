#!/bin/bash

# USAGE: bin/compare_md5.sh input_file.txt "expected_md5"
# this script will calculate the the md5sum for a file
# then compare it to the expected value
# if they do not match, the script exits with error code 1
# otherwise, do nothing (exits without error)


# check for helper script
if [ -f "$(dirname $0)/get_md5sum.sh" ]; then
    :
else
    echo "its not there" && exit 1
fi

# script args
input="${1:-none}"
[ "$input" == 'none' ] && printf 'input file not provided\n' && exit 1
[ ! -f "$input" ] && printf 'input item is not a file: ${input}\n' && exit 1
expected_md5="${2:-none}"
[ "$expected_md5" == 'none' ] && printf 'expected md5sum not provided\n' && exit 1

# get md5 for the input file
input_md5="$($(dirname $0)/get_md5sum.sh "${input}")"

# compare with the expected md5
if [ "${input_md5}" == "${expected_md5}" ]; then
    :
else
    echo "ERROR: the input md5 $input_md5 does not match the provided expected md5 $expected_md5"
    exit 1
fi
