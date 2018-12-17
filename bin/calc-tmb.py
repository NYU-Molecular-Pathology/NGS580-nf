#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calcalate the tumor mutation burden in variants/Megabase
"""
import sys
args = sys.argv
num_variants = float(args[1])
num_bases = float(args[2])
print( num_variants / num_bases * 1000000 )
