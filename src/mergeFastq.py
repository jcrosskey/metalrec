#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Merge 2 fastq files with paired end reads into one interleaved fastq file

Created on Fri Jan 30 10:41:02 EST 2015
@author: cjg
"""
import sys, os
import argparse
import time
import re
import itertools
## =================================================================
## Function: getRead
## =================================================================
def getReadFromFasta(fin):
    output = ""
    for line in fin:
        if line[0] != ">":
            output += line
        else:
            yield output
            output = line
    yield output


def getReadFromFastq(fin):
    output = ""
    count = 0
    for line in fin:
        if count != 4:
            output += line
            count = count + 1
        else:
            yield output
            output = line
            count = 1
    yield output


## =================================================================
## Function: splitReads
## =================================================================
def mergeFastq(Fastq1, Fastq2, out):
    ''' Merge 2 fastq files with paired end reads into one interleaved fastq file
        Input:  2 fastq files with paired end reads
        Output: interleaved merged file, default: sys.stdout
    '''
    count = 0
        
    for (r1,r2) in itertools.izip(getReadFromFastq(Fastq1), getReadFromFastq(Fastq2)):
        out.write('{}{}'.format(r1,r2))
        count = count + 1
        if count > 0 and count % 10000 == 0:
            sys.stderr.write("Merged {} reads..\n".formate(count))
    return count


## =================================================================
## argument parser
## =================================================================
parser = argparse.ArgumentParser(description="Split big equence file in fasta/fastq format into small ones",
                                 prog = 'splitReads', #program name
                                 prefix_chars='-', # prefix for options
                                 fromfile_prefix_chars='@', # if options are read from file, '@args.txt'
                                 conflict_handler='resolve', # for handling conflict options
                                 add_help=True, # include help in the options
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter # print default values for options in help message
                                 )

## input files and directories
parser.add_argument("fastqs", help = "2 fastq files to merge", nargs=2, type=argparse.FileType('r'))

## output directory
parser.add_argument("-o","--out",help="output fastq file",dest='outputFile',default=sys.stdout, type=argparse.FileType('w'))

parser.add_argument("-v","--verbose",help="verbose, more output",action='store_true',dest='verbose')
## =================================================================
## main function
## =================================================================
def main(argv=None):
    
    if argv is None:
        args = parser.parse_args()

    sys.stderr.write("\n===========================================================\n")
    start_time = time.time()

    seq_count = mergeFastq(args.fastqs[0], args.fastqs[1], args.outputFile)
    if args.verbose:
        sys.stderr.write("Total number of reads: {}\n".format(seq_count))

    sys.stderr.write("total time :" + str(time.time() - start_time) +  " seconds")
    sys.stderr.write("\n===========================================================\nDone\n")
##==============================================================
## call from command line (instead of interactively)
##==============================================================

if __name__ == '__main__':
    sys.exit(main())
