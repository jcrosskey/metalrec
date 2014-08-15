#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Clean up sam file generated from mapping, get rid of reads that are very noisy
Created on Fri Jun 18 13:14:18 2014
@author: cjg
"""
import sys
import argparse
import metalrec_lib
import time

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
parser.add_argument("-i","--in",help="input sam file",dest='samFile',required=True)
parser.add_argument("-i1","--in1",help="input ref sequence file",dest='seqFile',required=True)

parser.add_argument("--maxSub",help="maximum stretch of substitution",dest='maxSub',default=3, type=int)
parser.add_argument("--maxIns",help="maximum stretch of insertion",dest='maxIns',default=3, type=int)
parser.add_argument("--maxDel",help="maximum stretch of deletion",dest='maxDel',default=3, type=int)
parser.add_argument("--subRate",help="maximum substitution rate allowed",dest='maxSubRate',default=0.020, type=float)
parser.add_argument("--insRate",help="maximum insertion rate allowed",dest='maxInsRate',default=0.20, type=float)
parser.add_argument("--delRate",help="maximum deletion rate allowed",dest='maxDelRate',default=0.20, type=float)

## output directory
parser.add_argument("-o","--out",help="output file",dest='outputFile',required=True)

## =================================================================
## main function
## =================================================================
def main(argv=None):
    
    if argv is None:
        args = parser.parse_args()
    print "\n==========================================================="
    start_time = time.time()
    rSeq = metalrec_lib.read_single_seq(args.seqFile)
    metalrec_lib.clean_samfile(args.samFile, args.outputFile, rSeq, args.maxSub, args.maxIns, args.maxDel, args.maxSubRate, args.maxInsRate, args.maxDelRate)
    print "total time :" + str(time.time() - start_time) +  "seconds"
    print "==========================================================="
#    ref_bps, ref_ins_dict, consensus_seq, cov_depths = metalrec_lib.get_consensus(args.samFile,rSeq)
#    print "\n",consensus_seq,"\n"
#    print cov_depths

##==============================================================
## call from command line (instead of interactively)
##==============================================================

if __name__ == '__main__':
    sys.exit(main())

