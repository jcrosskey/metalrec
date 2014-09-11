#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
align_tmp.py

Align reads to corrected PacBio sequences at all correction steps.

Created on Fri Aug 22 10:46:18 2014
@author: cjg
"""
import sys
import os
import argparse
import time
import re
import glob

sys.path.append("/chongle/shared/software/metalrec/src")
import metalrec_lib
## =================================================================
## argument parser
## =================================================================
parser = argparse.ArgumentParser(description="Align reads to corrected PacBio sequences at all correction steps.",
                                 prog = 'align_tmp', #program name
                                 prefix_chars='-', # prefix for options
                                 fromfile_prefix_chars='@', # if options are read from file, '@args.txt'
                                 conflict_handler='resolve', # for handling conflict options
                                 add_help=True, # include help in the options
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter # print default values for options in help message
                                 )

## input files and directories
parser.add_argument("-f","--fasta",help="fasta file including all the good reads recruited by the PacBio sequence",dest='fastaFile',required=True)
parser.add_argument("-i","--inputDir",help="Directory where temporary directories reside",dest='inputDir',required=True)
## options for job queueing
parser.add_argument("-q","--queue",help="queue to use for jobs",dest='queue',default='large')
parser.add_argument("-p","--threads",help="number of threads to use every align job",dest='threads',default=2,type=int) # number of threads to use
parser.add_argument("-m","--memory",help="memory to request",dest='memory',default='8g') # 
## output directory
#parser.add_argument("-d","--outdir",help="output directory",dest='outputDir',default='./output/')

## =================================================================
## main function
## =================================================================
def main(argv=None):
    
    if argv is None:
        args = parser.parse_args()

    if not os.path.exists(args.inputDir):
        sys.exit("input directory does not exist!\n")

    args.inputDir = os.path.abspath(args.inputDir)
    args.fastaFile = os.path.abspath(args.fastaFile)
    inputDir = args.inputDir if args.inputDir[-1] == '/' else args.inputDir + '/'
    PacBio_seqs = glob.glob(inputDir + 'round*/seq.fasta')

    for PacBio_seq in PacBio_seqs:
        seq_dir = os.path.dirname(PacBio_seq) # directory for output of this subread
        # write bbmap script for this sequence
        bbmap_name = seq_dir + '/bbmap.sh'
        bbmap = open(bbmap_name,'w')
        bbmap.write('#!/bin/bash\n\n#PBS -l walltime=2:00:00\n#PBS -l nodes=1:ppn={}\n#PBS -q {}\n#PBS -N bbmap\n#PBS -e {}\n#PBS -o {}\n'.format(str(args.threads), args.queue, seq_dir+'/bbmap.err',seq_dir+'/bbmap.out' ))
        bbmap.write("\ncd {}\n".format(seq_dir))
        bbmap.write("echo Starting Time is $(date)\n")
        bbmap.write("bbmapskimmer.sh build=1 ref={}\n\n".format(PacBio_seq))
        bbmap.write("bbmapskimmer.sh outputunmapped=f build=1 killbadpairs=f ambiguous=all local=f maxindel=5 maxindel2=50 strictmaxindel=t maxsublen=3 keepnames=f  minid=0.70 k=10 ignorebadquality=t secondary=t maxsites=50 sam=1.4 requirecorrectstrand=f idtag=t saa=f md=t {} threads={} trimreaddescriptions=t in={} out=bbmap.sam \n\n".format('-Xmx'+args.memory, str(args.threads), seq_dir + '/reads.fasta'))
        bbmap.write("bbmapskimmer.sh outputunmapped=f build=1 killbadpairs=f ambiguous=all local=f maxindel=5 maxindel2=50 strictmaxindel=t maxsublen=3 keepnames=f  minid=0.70 k=10 ignorebadquality=t secondary=t maxsites=50 sam=1.4 requirecorrectstrand=f idtag=t saa=f md=t {} threads={} trimreaddescriptions=t in={} out=bbmap_allreads.sam \n\n".format('-Xmx'+args.memory, str(args.threads), args.fastaFile))
        #bbmap.write("python /chongle/shared/software/metalrec/src/clean_sam.py -i bbmap.sam -i1 seq.fasta --maxSub 100 --maxIns 100 --maxDel 100 --subRate 0.1 --insRate 0.2 --delRate 0.2 -o realign.sam\n")
        bbmap.write("echo Ending Time is $(date)\n")
        bbmap.close()
        os.system("qsub {}".format(bbmap_name)) #run system command and submit bbmap job

##==============================================================
## call from command line (instead of interactively)
##==============================================================

if __name__ == '__main__':
    sys.exit(main())

