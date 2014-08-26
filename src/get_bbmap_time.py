#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
get_bbmap_time.py

Get the bbmap running time for the mapping jobs.
Created on Tue Aug 26 11:05:03 EDT 2014
@author: cjg
"""
import sys
import os
import argparse
import time
import re
import glob

#sys.path.append("/chongle/shared/software/metalrec/src")
#import metalrec_lib

## =================================================================
## argument parser
## =================================================================
parser = argparse.ArgumentParser(description="Get the bbmap running time for the mapping jobs.",
                                 prog = 'check_PBreads', #program name
                                 prefix_chars='-', # prefix for options
                                 fromfile_prefix_chars='@', # if options are read from file, '@args.txt'
                                 conflict_handler='resolve', # for handling conflict options
                                 add_help=True, # include help in the options
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter # print default values for options in help message
                                 )

## input files and directories
parser.add_argument("-i","--inputDir",help="input directory",dest='inputDir',required=True)
## output files and directories
parser.add_argument("-o","--out",help="output file",dest='outputFile',required=True)

## =================================================================
def get_time(inputfile):
    if not os.path.exists(inputfile):
        sys.exit('{} does not exist!\n'.format(inputfile))
    times = []
    with open(inputfile,'r') as fin:
        for line in fin:
            if 'Total time' in line:
                times.append(re.findall("\d+.\d+",line)[0])
    if len(times) != 2:
        sys.stdout.write("Number of lines with Total time is {}, not 2!\n".format(len(times)))
    return times
## =================================================================
## main function
## =================================================================
def main(argv=None):
    if argv is None:
        args = parser.parse_args()

    if not os.path.exists(os.path.dirname(os.path.abspath(args.outputFile))):
        os.makedirs(os.path.dirname(os.path.abspath(args.outputFile)))
    fout = open(args.outputFile,'w')
    fout.write('read_name\tindex_time\tmapping_time\n')
    err_files = glob.glob(os.path.abspath(args.inputDir) + '/*/bbmap.err')
    sys.stdout.write('Number of err files to read: {}\n'.format(len(err_files)))
    for err_file in err_files:
        read_name = os.path.abspath(err_file).split('/')[-2]
        times = get_time(err_file)
        if len(times) == 2:
            fout.write('{}\t{}\t{}\n'.format(read_name, times[0], times[1]))
        else:
            fout.write('{}\t{}\t{}\n'.format(read_name, 'NA', 'NA'))
    fout.close()
    
##==============================================================
## call from command line (instead of interactively)
##==============================================================

if __name__ == '__main__':
    sys.exit(main())

