#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
test drive
Created on Fri May 23 13:14:18 2014
@author: cjg
"""
import sys
import argparse
import metalrec_lib

## =================================================================
## argument parser
## =================================================================
parser = argparse.ArgumentParser(description="parse sam file and get summary statistics",
                                 prog = 'testDrive', #program name
                                 prefix_chars='-', # prefix for options
                                 fromfile_prefix_chars='@', # if options are read from file, '@args.txt'
                                 conflict_handler='resolve', # for handling conflict options
                                 add_help=True, # include help in the options
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter # print default values for options in help message
                                 )

## input files and directories
parser.add_argument("-i","--in",help="input file",dest='inputFile',required=True)

## output directory
parser.add_argument("-o","--out",help="output file",dest='outputFile',required=True)

## =================================================================
## main function
## =================================================================
def main(argv=None):
    
    if argv is None:
        args = parser.parse_args()

    metalrec_lib.rm_bad_record(args.inputFile, args.outputFile)

##==============================================================
## call from command line (instead of interactively)
##==============================================================

if __name__ == '__main__':
    sys.exit(main())
