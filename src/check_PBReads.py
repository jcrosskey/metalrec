#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Check the PacBio subreads and see if any of them has good regions
Created on Thu Jul 04 14:33:18 2014
@author: cjg
"""
import sys
import os
import argparse
import time
import re

sys.path.append("/chongle/shared/software/metalrec/src")
import metalrec_lib
## =================================================================
## argument parser
## =================================================================
parser = argparse.ArgumentParser(description="Check the PacBio subreads and see if any of them has good regions",
                                 prog = 'check_PBreads', #program name
                                 prefix_chars='-', # prefix for options
                                 fromfile_prefix_chars='@', # if options are read from file, '@args.txt'
                                 conflict_handler='resolve', # for handling conflict options
                                 add_help=True, # include help in the options
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter # print default values for options in help message
                                 )

## input files and directories
parser.add_argument("-f","--fasta",help="input PacBio filtered subreads sequence file",dest='fastaFile',required=True)
parser.add_argument("-pe",help="paired end Illumina reads data set",dest='peFiles', nargs='+') # paired end reads files, in pair of separate files
parser.add_argument("-se",help="paired end Illumina reads data set",dest='seFiles', nargs='+') # single end reads files. TODO: what about interleaved PE reads?
## output directory
parser.add_argument("-d","--outdir",help="output directory",dest='outputDir',default='./output/')

## options to control how many sequences to process
parser.add_argument("-c","--count",help="number of subreads to process",dest='count',default=-1, type=int)
parser.add_argument("-k","--skip",help="number of subreads to skip",dest='skip',default=0, type=int)
parser.add_argument("-e","--last",help="index of the last subread to look at",dest='last',default=0, type=int)

## options for job queueing
parser.add_argument("-q","--queue",help="queue to use for jobs",dest='queue',default='large') # queue to use
parser.add_argument("-p","--threads",help="number of threads to use every align job",dest='threads',default=2,type=int) # number of threads to use
parser.add_argument("-m","--memory",help="memory to request",dest='memory',default='8g') # memory to use

## setting thresholds
#parser.add_argument("--maxSub",help="maximum stretch of substitution",dest='maxSub',default=3, type=int)
#parser.add_argument("--maxIns",help="maximum stretch of insertion",dest='maxIns',default=3, type=int)
#parser.add_argument("--maxDel",help="maximum stretch of deletion",dest='maxDel',default=3, type=int)
#parser.add_argument("--subRate",help="maximum substitution rate allowed",dest='maxSubRate',default=0.020, type=float)
#parser.add_argument("--insRate",help="maximum insertion rate allowed",dest='maxInsRate',default=0.20, type=float)
#parser.add_argument("--delRate",help="maximum deletion rate allowed",dest='maxDelRate',default=0.20, type=float)
#parser.add_argument("--minCV",help="minimum coverage depth",dest='minCV',default=1, type=int)
parser.add_argument("--minPacBioLen",help="minimum length of PacBio read to be aligned",dest='minPacBioLen',default=1000, type=int)
#parser.add_argument("--minGoodLen",help="minimum contiguous well covered region length",dest='minGoodLen',default=1000, type=int)
#parser.add_argument("--polyN",help="minimum number of read support for a base to be considered",dest='minReads',default=3, type=int)
#parser.add_argument("--polyR",help="minimum proportion of read support for a base to be considered",dest='minPercent',default=0.01, type=float)

#parser.add_argument("-o","--out",help="output file",dest='outputFile',required=True)

## =================================================================
## main function
## =================================================================
def main(argv=None):
    
    if argv is None:
        args = parser.parse_args()

    Illumina_in = '' # 'in=...' string part in bbmap.sh command
    if args.peFiles is not None: # if there are PE reads data
        in1 = ''
        in2 = ''
        for i in xrange(len(args.peFiles)/2):
            if i != 0:
                in1 += args.peFiles[i*2]
                in2 += args.peFiles[i*2 + 1]
            else:
                in1 += ',' + args.peFiles[i*2]
                in2 += ',' + args.peFiles[i*2 + 1]
        if args.seFiles is not None:
            for i in xrange(len(args.seFiles)):
                in1 += ',' + args.seFiles[i]
                in2 += ',null'
        Illumina_in = 'in1=' + in1 + ' in2= ' + in2 # concatenate in1= and in2= together
        map_out = 'out=' + ','.join(['bbmap.sam'] * (len(args.peFiles)/2 + len(args.seFiles))) # out= string in bbmap.sh
    else: # if there are only SE reads data
        if args.seFiles is None:
            sys.exit("no input Illumina reads specified!")
        else:
            Illumina_in = ','.join(args.seFiles)
            Illumina_in = 'in=' + Illumina_in
        map_out = 'out=' + ','.join(['bbmap.sam'] *  len(args.seFiles))

    # make sure output directory exists
    if not os.path.exists(args.outputDir):
        os.makedirs(args.outputDir)
    args.outputDir = os.path.abspath(args.outputDir) # make sure the path is the full path

    read_count = 0 # number of reads that are processed (not including skipped reads)
    scan_count = 0 # number of reads that are scanned
    skip_count = 0 # number of reads that are skipped
    read_seq = ''
    skip = False

    # scan the fasta file for PacBio subreads
    fsr = open(args.fastaFile,'r')
    line = fsr.readline()
    nextline = ''
    while line != '' and (args.count == -1 or read_count < args.count) and (scan_count < args.last or args.last == 0):
        if line[0] == '>':
            read_name = line.split()[0][1:] # get read name
            read_seq = ''
            scan_count += 1 # increase scanned read count

            # check if enough reads are skipped
            if skip_count < args.skip : # skip certain number of reads
                skip = True
                skip_count += 1
                #sys.stdout.write("{}\t{}\t{}\t{}\tskip\n".format(scan_count,skip_count, read_count,read_name))
                nextline = fsr.readline().strip('\n')
                while nextline[0] != '>':
                    nextline = fsr.readline().strip('\n')
                line = nextline
                continue
            else: # first read ('') or skipped enough reads
                skip = False
            
            if not skip:
                sys.stdout.write("{}\t{}\t{}\t{}".format(scan_count,skip_count, read_count,read_name))
                legal_name = re.sub('/','__',read_name)
                seq_dir = args.outputDir + '/' + legal_name # directory for output of this subread
                # first check if this read is already done
                if not os.path.exists(seq_dir+'/bbmap.sam'): # if this read is not processed yet. criterion for already processed: bbmap.sam exists
                    read_count += 1 # increase processed read count 
                    nextline = fsr.readline().strip('\n')
                    while nextline != '' and nextline[0] != '>':
                        read_seq += nextline
                        nextline = fsr.readline().strip('\n')
                    line = nextline
                    sys.stdout.write("\t{}".format(len(read_seq)))
                    if len(read_seq) >= args.minPacBioLen : # if sequence length passes length threshold, write the fasta file
                        if not os.path.exists(seq_dir):
                            os.makedirs(seq_dir)
                        sys.stdout.write("\t processing\n")
                        # write fasta file for this sequence
                        fasta_name = seq_dir + '/' + legal_name + '.fasta'
                        myfasta = open(fasta_name,'w')
                        myfasta.write('>{}\n{}'.format(read_name,read_seq))
                        myfasta.close()

                        # write bbmap script for this sequence
                        bbmap_name = seq_dir + '/bbmap.sh'
                        bbmap = open(bbmap_name,'w')
                        bbmap.write('#!/bin/bash\n\n#PBS -l walltime=10:00:00\n#PBS -l nodes=1:ppn={}\n#PBS -q {}\n#PBS -N bbmap\n#PBS -e {}\n#PBS -o {}\n'.format(str(args.threads), args.queue, seq_dir+'/bbmap.err',seq_dir+'/bbmap.out' ))
                        bbmap.write("\ncd {}\n".format(seq_dir))
                        bbmap.write("echo Starting Time is $(date)\n")
                        bbmap.write("bbmapskimmer.sh build=1 ref={} k=10 {}\n\n".format(fasta_name, '-Xmx'+args.memory))
                        bbmap.write("bbwrap.sh mapper=bbmappacbioskimmer  outputunmapped=f build=1 killbadpairs=f ambiguous=all local=f strictmaxindel=t keepnames=t k=10 ignorebadquality=t secondary=t maxsites=50 sam=1.4 requirecorrectstrand=f idtag=t saa=f md=t {} threads={} trimreaddescriptions=t {} {} append\n".format('-Xmx'+args.memory, str(args.threads), Illumina_in, map_out))
                        bbmap.write("echo Ending Time is $(date)\n")
                        bbmap.close()
                        time.sleep(0.5)
                        os.system("qsub {}".format(bbmap_name)) #run system command and submit bbmap job
                    else:
                        sys.stdout.write("\t too short\n")
                else:
                    sys.stdout.write("\t already done\n")
                    nextline = fsr.readline().strip('\n')
                    while nextline[0] != '>':
                        nextline = fsr.readline().strip('\n')
                    line = nextline
                    #read_seq = ''
                    #read_name = line.split()[0][1:] # filtered subread name
                    #continue

##==============================================================
## call from command line (instead of interactively)
##==============================================================

if __name__ == '__main__':
    sys.exit(main())

