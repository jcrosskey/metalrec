#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
run_metalrec

Created on Wed, Aug 20, 14:06:18 2014
@author: chaij@ornl.gov
"""
import sys, os
import argparse
import metalrec_lib
import time
import glob
import re

## =================================================================
## argument parser
## =================================================================
parser = argparse.ArgumentParser(description="Error correction of PacBio reads from aligned results",
                                 prog = 'run_metalrec', #program name
                                 prefix_chars='-', # prefix for options
                                 fromfile_prefix_chars='@', # if options are read from file, '@args.txt'
                                 conflict_handler='resolve', # for handling conflict options
                                 add_help=True, # include help in the options
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter # print default values for options in help message
                                 )

## input files and directories
parser.add_argument("-i","--in",help="input directory",dest='inputDir',required=True)

## output directory
#parser.add_argument("-o","--out",help="output corrected PacBio sequence file",dest='oSeqFile',default=None)

# options
#parser.add_argument("-m",help="read mode, s(ingle) or p(air)",dest='rMode',default='p',choices=['s','p'])
#parser.add_argument("-v","--verbose",help="verbose, more output",action='store_true',dest='verbose')
parser.add_argument("-q","--queue",help="queue to use for jobs",dest='queue',default='large') # queue to use
parser.add_argument("-t","--wtime",help="wall time",dest='wtime',default='00:20:00') # wall time

## setting thresholds
parser.add_argument("--maxSub",help="maximum stretch of substitution",dest='maxSub',default=100, type=int)
parser.add_argument("--maxIns",help="maximum stretch of insertion",dest='maxIns',default=100, type=int)
parser.add_argument("--maxDel",help="maximum stretch of deletion",dest='maxDel',default=100, type=int)
parser.add_argument("--subRate",help="maximum substitution rate allowed",dest='maxSubRate',default=0.10, type=float)
parser.add_argument("--insRate",help="maximum insertion rate allowed",dest='maxInsRate',default=0.20, type=float)
parser.add_argument("--delRate",help="maximum deletion rate allowed",dest='maxDelRate',default=0.20, type=float)
parser.add_argument("--minCV",help="minimum coverage depth",dest='minCV',default=1, type=int)
parser.add_argument("--minPacBioLen",help="minimum PacBio length to be considered",dest='minPacBioLen',default=1000, type=int)
parser.add_argument("--minGoodLen",help="minimum contiguous well covered region length",dest='minGoodLen',default=400, type=int)
parser.add_argument("--polyN",help="minimum number of read support for a base to be considered",dest='minReads',default=3, type=int)
parser.add_argument("--polyR",help="minimum proportion of read support for a base to be considered",dest='minPercent',default=0.3, type=float)
## =================================================================
## main function
## =================================================================
def main(argv=None):
    # parse command line arguments    
    if argv is None:
        args = parser.parse_args()

    # check input and output file settings
    if not os.path.exists(args.inputDir):
        sys.exit("input directory {} does not exist!\n".format(args.inputDir))
    if args.inputDir[-1] != '/':
        args.inputDir = args.inputDir + '/'

    seqDirs = glob.glob(args.inputDir + '*')
    for seqDir in seqDirs: # for every PacBio sequence (directory with alignments)
        seqDir = os.path.abspath(seqDir)
        seqFile = glob.glob(seqDir + '/*.fasta')[0]
        samFile = glob.glob(seqDir + '/bbmap.sam')[0]
        # write pbcorrect script for this sequence
        pbcorrect_name = seqDir + '/pbcorrect.sh'
        pbcorrect_out = seqDir + '/pbcorrect.out'
        #if os.path.exists(pbcorrect_name) and os.path.exists(pbcorrect_out) and os.path.exists(seqDir + '/corrected_PacBio.fasta'):
        if os.path.exists(pbcorrect_name) and os.path.exists(pbcorrect_out):
            sys.stdout.write("{} already done\n".format(seqDir))
        else:
            pbcorrect = open(pbcorrect_name,'w')
            pbcorrect.write('#!/bin/bash\n\n#PBS -l walltime={}\n#PBS -l nodes=1:ppn={}\n#PBS -q {}\n#PBS -N pbcorrect\n#PBS -e {}\n#PBS -o {}\n'.format(args.wtime, 1, args.queue, seqDir+'/pbcorrect.err',seqDir+'/pbcorrect.out' ))
            pbcorrect.write("\ncd {}\n".format(seqDir))
            pbcorrect.write("echo Starting Time is $(date)\n")
            pbcorrect.write("python /chongle/shared/software/metalrec/src/metalrec.py  -i {} -s {} --maxSub {} --maxIns {} --maxDel {} --subRate {} --insRate {} --delRate {}  --minCV {} --minPacBioLen {} --minGoodLen {} --polyN {} --polyR {}\n\n".format(seqFile, samFile, args.maxSub, args.maxIns, args.maxDel, args.maxSubRate, args.maxInsRate, args.maxDelRate, args.minCV, args.minPacBioLen, args.minGoodLen, args.minReads, args.minPercent))
            pbcorrect.write("echo Ending Time is $(date)\n")
            pbcorrect.close()
            time.sleep(0.1)
            os.system("qsub {}".format(pbcorrect_name)) #run system command and submit pbcorrect job
##==============================================================
## call from command line (instead of interactively)
##==============================================================

if __name__ == '__main__':
    sys.exit(main())
