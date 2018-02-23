#!/bin/bash

# USAGE: get_md5sum.sh input_file
# OUTPUT: md5sum for the file

# get the md5 hash for a file on Linux or macOS / OS X
# script to checks whether 'md5sum' or 'md5' programs are installed,
# e.g. if running on macOS / OS X or Linux, etc
# and runs the appropriate commands
# and prints md5 hash to console

# $ md5sum --version
# md5sum (GNU coreutils) 8.4
# Copyright (C) 2010 Free Software Foundation, Inc.
# License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.
# This is free software: you are free to change and redistribute it.
# There is NO WARRANTY, to the extent permitted by law.
# Written by Ulrich Drepper, Scott Miller, and David Madore.

# using md5 version from macOS 10.12.5

# script args
input="${1:-none}"
[ "$input" == 'none' ] && printf 'at least one argument is required\n' && exit 1

# check for programs & get hash
if type md5sum > /dev/null 2>&1; then
    md5_sum="$(md5sum "$input" | cut -d ' ' -f1)"

elif type md5 > /dev/null 2>&1; then
    md5_sum="$(md5 -q "$input")"

else
    printf 'neither md5 nor md5sum are installed\n'
    exit 1
fi

printf '%s' "$md5_sum"
exit
