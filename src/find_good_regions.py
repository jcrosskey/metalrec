#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Find good regions for a PacBio read
Created on Thu Jul 03 13:49:18 2014
@author: cjg
"""
import sys
import argparse
import metalrec_lib
import time

## =================================================================
## argument parser
## =================================================================
parser = argparse.ArgumentParser(description="Find good regions of a PacBio read using mapping from short reads onto it",
                                 prog = 'find_good_regions', #program name
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
parser.add_argument("--minCV",help="minimum coverage depth",dest='minCV',default=10, type=int)
parser.add_argument("--minPacBioLen",help="minimum contiguous well covered region length",dest='minPacBioLen',default=1000, type=int)
parser.add_argument("--polyN",help="minimum number of read support for a base to be considered",dest='minReads',default=3, type=int)
parser.add_argument("--polyR",help="minimum proportion of read support for a base to be considered",dest='minPercent',default=0.01, type=float)


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
    rname = 'unknown'
    with open(args.samFile, 'r') as mysam:
        for line in mysam:
            if line[0] == '@': # header line
                if line[1:3] == 'SQ': # reference sequence dictionary
                    rname = line[(line.find('SN:') + len('SN:')) : line.find('\t',line.find('SN:'))] # reference sequence name
                    break
    ref_bps, ref_ins_dict, readinfo = metalrec_lib.read_and_process_sam(args.samFile, rSeq, args.maxSub, args.maxIns, args.maxDel, args.maxSubRate, args.maxInsRate, args.maxDelRate,args.minPacBioLen, args.minCV)
    good_regions = metalrec_lib.get_good_regions(ref_bps, rSeq, args.minPacBioLen, args.minCV)
    sys.stdout.write("good regions" + str(good_regions) + "\n")

    myout = open(args.outputFile,'a')
    myout.write(rname)
    if len(good_regions) > 0:
        for region in good_regions:
            myout.write("\t{}".format(str(region)))
    myout.write('\n')
    myout.close()
    print "total time :" + str(time.time() - start_time) +  " seconds"
    print "==========================================================="

##==============================================================
## call from command line (instead of interactively)
##==============================================================

if __name__ == '__main__':
    sys.exit(main())

