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
parser.add_argument("-c","--count",help="number of subreads to process",dest='count',default=-1, type=int)
parser.add_argument("-k","--skip",help="number of subreads to skip",dest='skip',default=0, type=int)

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
parser.add_argument("-d","--outdir",help="output directory",dest='outputDir',default='./output/')
#parser.add_argument("-o","--out",help="output file",dest='outputFile',required=True)

## =================================================================
## main function
## =================================================================
def main(argv=None):
    
    if argv is None:
        args = parser.parse_args()

    # Illumina data
    pe800="/chongle/shared/database/03_PacBio/Wetlands/Illumina/800bp_pe.fastq"
    se800="/chongle/shared/database/03_PacBio/Wetlands/Illumina/800bp_se.fastq"
    pe400="/chongle/shared/database/03_PacBio/Wetlands/Illumina/400bp_pe.fastq"
    se400="/chongle/shared/database/03_PacBio/Wetlands/Illumina/400bp_se.fastq"
    pe270="/chongle/shared/database/03_PacBio/Wetlands/Illumina/270bp_pe.fastq"
    se270="/chongle/shared/database/03_PacBio/Wetlands/Illumina/270bp_se.fastq"

    # make sure output directory exists
    if not os.path.exists(args.outputDir):
        os.makedirs(args.outputDir)

    read_count = 0 # number of reads that are processed (not including skipped reads)
    scan_count = 0 # number of reads that are scanned
    skip_count = 0 # number of reads that are skipped
    read_seq = ''
    skip = False

    # scan the fasta file for PacBio subreads
    fsr = open(args.fastaFile,'r')
    line = fsr.readline()
    nextline = ''
    while line != '' and (args.count == -1 or read_count < args.count):
        if line[0] == '>':
            read_name = line.split()[0][1:] # get read name
            read_seq = ''
            scan_count += 1 # increase scanned read count

            # check if enough reads are skipped
            if skip_count < args.skip : # skip certain number of reads
                skip = True
                skip_count += 1
                sys.stdout.write("{}\t{}\t{}\t{}\tskip\n".format(scan_count,skip_count, read_count,read_name))
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
                if not os.path.exists(seq_dir): # if this read is not processed yet
                    os.makedirs(seq_dir)
                    read_count += 1 # increase processed read count 
                    nextline = fsr.readline().strip('\n')
                    while nextline != '' and nextline[0] != '>':
                        read_seq += nextline
                        nextline = fsr.readline().strip('\n')
                    line = nextline
                    sys.stdout.write("\t{}".format(len(read_seq)))
                    if len(read_seq) >= args.minPacBioLen : # if sequence length passes length threshold, write the fasta file
                        sys.stdout.write("\t processing\n")
                        # write fasta file for this sequence
                        fasta_name = seq_dir + '/' + legal_name + '.fasta'
                        myfasta = open(fasta_name,'w')
                        myfasta.write('>{}\n{}'.format(read_name,read_seq))
                        myfasta.close()

                        # write bbmap script for this sequence
                        bbmap_name = seq_dir + '/bbmap.sh'
                        bbmap = open(bbmap_name,'w')
                        bbmap.write('#!/bin/bash\n\n#PBS -l walltime=10:00:00\n#PBS -l nodes=1:ppn=16\n#PBS -q large\n#PBS -N bbmap\n#PBS -e {}\n#PBS -o {}\n'.format(seq_dir+'/bbmap.err',seq_dir+'/bbmap.out' ))
                        bbmap.write("\ncd {}\n".format(seq_dir))
                        bbmap.write("echo Starting Time is $(date)\n")
                        bbmap.write("bbmapskimmer.sh build=1 ref={}\n".format(fasta_name))
                        bbmap.write("bbwrap.sh mapper=bbmappacbioskimmer  outputunmapped=f build=1 killbadpairs=f ambiguous=all local=f maxindel=5 maxindel2=50 strictmaxindel=t maxsublen=3 keepnames=t  minid=0.70 k=10 ignorebadquality=t secondary=t maxsites=50 sam=1.4 requirecorrectstrand=f idtag=t saa=f md=t -Xmx50g threads=16 trimreaddescriptions=t in={},{} out=bbmap.sam,bbmap.sam append\n".format(pe800, se800))
                        bbmap.write("echo Ending Time is $(date)\n")
                        bbmap.close()
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

