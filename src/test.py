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
parser.add_argument("-i","--in",help="input sam file",dest='samFile',required=True)
parser.add_argument("-i1","--in1",help="input ref sequence file",dest='seqFile',required=True)

## output directory
parser.add_argument("-o","--out",help="output file",dest='outputFile',required=True)
parser.add_argument("-o1","--out1",help="single end output file",dest='soutputFile',required=True)

## =================================================================
## main function
## =================================================================
def main(argv=None):
    
    if argv is None:
        args = parser.parse_args()

    rSeq = metalrec_lib.read_single_seq(args.seqFile)
    ref_bps, ref_ins_dict, readinfo = metalrec_lib.read_and_process_sam_pair(args.samFile,rSeq,maxSubRate=0.2, outsam = args.outputFile)
    print len(readinfo)
    ref_bps, ref_ins_dict, readinfo = metalrec_lib.read_and_process_sam(args.samFile,rSeq,maxSubRate=0.2, outsam = args.soutputFile)
    print len(readinfo)
    #ref_bps, ref_ins_dict, consensus_seq, cov_depths = metalrec_lib.get_consensus(args.samFile,rSeq)
    #print "\n",consensus_seq,"\n"
    #print cov_depths
    ##metalrec_lib.rm_bad_record(args.samFile,args.outputFile)

##==============================================================
## call from command line (instead of interactively)
##==============================================================

if __name__ == '__main__':
    sys.exit(main())
