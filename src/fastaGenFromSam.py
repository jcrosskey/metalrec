#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
fastaGenFromSam

Created on Tue Oct 28 09:10:01 EDT 2014

@author: cjg
"""
import sys
import argparse
import samread 

## =================================================================
## generate fasta file including mapped sequences from sam file
## =================================================================
def fastaGenFromSam(samFile, out, checkReads = True, verbose=False):
    ''' generate fasta file including mapped sequences from sam file
        Input:  samFile - sam file from mapping/alignment result
        Output: out - file object (where the sequences will be written to)
    '''
    recorded_reads = []
    mapping_count = 0
    with open(samFile, 'r') as sam:
        for line in sam:
            if line[0] == '@':
                continue
            else:
                myread = samread.SamRead(line)
                if checkReads and myread.is_unmapped():
                    continue
                mapping_count += 1
                if myread.qname not in recorded_reads:
                    out.write(">{}\n{}\n".format(myread.qname, myread.qSeq))
                    recorded_reads.append(myread.qname)
    if verbose:
        sys.stdout.write("{} reads with {} mappings\n".format(len(recorded_reads), mapping_count))
    

## =================================================================
## argument parser
## =================================================================
parser = argparse.ArgumentParser(description="generate fasta file including mapped sequences from sam file",
                                 prog = 'fastaGenFromSam', #program name
                                 prefix_chars='-', # prefix for options
                                 fromfile_prefix_chars='@', # if options are read from file, '@args.txt'
                                 conflict_handler='resolve', # for handling conflict options
                                 add_help=True, # include help in the options
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter # print default values for options in help message
                                 )

## input files and directories
parser.add_argument("-i","--in",help="input sam file",dest='samFile',required=True)

## output directory
parser.add_argument("-o","--out",help="output fasta file",dest='outputFile',default=sys.stdout, type = argparse.FileType('w'))

## options
parser.add_argument("-v", "--verbose", help="verbose mode", dest = "verbose", action="store_true")
parser.add_argument("-c", "--checkReads", help="check and see if read is mapped before writing into fasta", dest = "checkReads", action="store_true")

## =================================================================
## main function
## =================================================================
def main(argv=None):
    if argv is None:
        args = parser.parse_args()

    fastaGenFromSam(args.samFile, args.outputFile, args.checkReads, args.verbose)

##==============================================================
## call from command line (instead of interactively)
##==============================================================

if __name__ == '__main__':
    sys.exit(main())
