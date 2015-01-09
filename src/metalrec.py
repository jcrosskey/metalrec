#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
metalrec

Created on Wed, Aug 20, 14:06:18 2014
@author: chaij@ornl.gov
"""
import sys, os
import argparse
import metalrec_lib
import numpy
import time
import re

## =================================================================
## argument parser
## =================================================================
parser = argparse.ArgumentParser(description="Error correction of a PacBio read from alignment of Illumina reads to this PacBio read.",
                                 prog = 'metalrec', #program name
                                 prefix_chars='-', # prefix for options
                                 fromfile_prefix_chars='@', # if options are read from file, '@args.txt'
                                 conflict_handler='resolve', # for handling conflict options
                                 add_help=True, # include help in the options
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter # print default values for options in help message
                                 )

## input files and directories
parser.add_argument("-i","--in",help="input ref sequence file",dest='seqFile',required=True)
parser.add_argument("-s","--sam",help="input sam file",dest='samFile',required=True)

## output directory
parser.add_argument("-o","--out",help="output corrected PacBio sequence file",dest='oSeqFile',default=None)
parser.add_argument("-od","--outDir",help="directory for the intermediate files",dest='outDir',default = None)

# options
#parser.add_argument("-m",help="read mode, s(ingle) or p(air)",dest='rMode',default='p',choices=['s','p'])
parser.add_argument("-v","--verbose",help="verbose, more output",action='store_true',dest='verbose')
parser.add_argument("--width",help="print width for sequences in verbose output",dest='width', default = 100, type = int)
parser.add_argument("--checkEnds",help="check the substitution errors at the ends",dest='checkEnds', action='store_false')

## setting thresholds
parser.add_argument("--minOverlap",help="minimum overlap length between reads",dest='minOverlap',default=10, type=int)
parser.add_argument("--minOverlapRatio",help="minimum ratio of average overlap length between reads",dest='minOverlapRatio',default=0.1, type=float)
parser.add_argument("--maxSub",help="maximum stretch of substitution",dest='maxSub',default=-1, type=int)
parser.add_argument("--maxIns",help="maximum stretch of insertion",dest='maxIns',default=-1, type=int)
parser.add_argument("--maxDel",help="maximum stretch of deletion",dest='maxDel',default=-1, type=int)
parser.add_argument("--subRate",help="maximum substitution rate allowed",dest='maxSubRate',default=0.05, type=float)
parser.add_argument("--indelRate",help="maximum insertion rate allowed",dest='maxInDelRate',default=0.30, type=float)
parser.add_argument("--minCV",help="minimum coverage depth",dest='minCV',default=1, type=int)
parser.add_argument("--minPacBioLen",help="minimum PacBio length to be considered",dest='minPacBioLen',default=1000, type=int)
parser.add_argument("--minGoodLen",help="minimum contiguous well covered region length",dest='minGoodLen',default=400, type=int)
## =================================================================
## main function
## =================================================================
def main(argv=None):
    # parse command line arguments    
    if argv is None:
        args = parser.parse_args()

    sys.stderr.write("\n===========================================================\n")
    start_time = time.time()
    # check input and output file settings
    ## required input files: PacBio sequence file and sam file for this PacBio sequence
    if not os.path.exists(args.seqFile):
        sys.exit("input PacBio sequence file does not exist!\n")
    if not os.path.exists(args.samFile):
        sys.exit("input sam file does not exist!\n")

    ## output file and directories, optional
    if args.oSeqFile is None: # default destination for the corrected PacBio sequence(contigs if the sequence is split into different regions)
        args.oSeqFile = os.path.dirname(os.path.abspath(args.seqFile))+ '/EC.fasta'
    #shortSeqFile = args.oSeqFile + '.short'
    if os.path.exists(args.oSeqFile): # overwrite the output file if it already exists
        os.remove(args.oSeqFile)
        sys.stderr.write("Output sequence file already exists, overwrite.\n")
    #if os.path.exists(shortSeqFile): # overwrite the output file if it already exists
    #    os.remove(shortSeqFile)
    elif not os.path.exists(os.path.dirname(os.path.abspath(args.oSeqFile))): # make sure the directory for the output file exists
        os.makedirs(os.path.dirname(os.path.abspath(args.oSeqFile)))

    if args.verbose:
        if args.outDir is None:
            args.outDir = os.path.dirname(os.path.abspath(args.samFile)) + '/EC/'
        else:
            args.outDir = os.path.abspath(args.outDir) + '/'
        if not os.path.exists(args.outDir):
            os.makedirs(args.outDir)
        sys.stderr.write("verbose output directory: {}.\n".format(args.outDir))

    if args.verbose:
        sys.stderr.write("minimum overlap length: {}\n".format(args.minOverlap))
        sys.stderr.write("minimum overlap length ratio: {}\n".format(args.minOverlapRatio))
        sys.stderr.write("maximum stretch of substitution: {}\n".format(args.maxSub))
        sys.stderr.write("maximum stretch of insertion: {}\n".format(args.maxIns))
        sys.stderr.write("maximum stretch of deletion: {}\n".format(args.maxDel))
        sys.stderr.write("maximum substitution rate allowed: {}\n".format(args.maxSubRate))
        sys.stderr.write("maximum indel rate allowed: {}\n".format(args.maxInDelRate))
        sys.stderr.write("minimum coverage depth: {}\n".format(args.minCV))
        sys.stderr.write("minimum PacBio read length to be considered: {}\n".format(args.minPacBioLen))
        sys.stderr.write("minimum good region length: {}\n".format(args.minGoodLen))
        sys.stderr.write("verbose mode: {}\n".format(args.verbose))
        sys.stderr.write("check substitution error rate at ends: {}\n".format(args.checkEnds))

    # read the PacBio sequence into memory
    rseq = metalrec_lib.read_single_seq(args.seqFile)
    # process sam file and save the read info
    s_time = time.time()
    ref_bps, ref_ins_dict, read_info = metalrec_lib.read_and_process_sam_samread(args.samFile, rseq, maxSub=args.maxSub, maxIns=args.maxIns, maxDel=args.maxDel,maxSubRate=args.maxSubRate, maxInDelRate=args.maxInDelRate, minPacBioLen=args.minPacBioLen, checkEnds=args.checkEnds, outDir=args.outDir, verbose=args.verbose)
    e_time = time.time()
    sys.stderr.write("processing sam file time :" + str(e_time - s_time) +  " second\n")

    if len(ref_bps) == 0: # empty sam file, or nothing
        sys.stderr.write("PacBio read does not have any coverage from Illumina reads\n")
    else:
        good_regions, cov_bps, avg_cov_depth = metalrec_lib.get_good_regions(ref_bps, rseq, minGoodLen=args.minGoodLen, minCV=args.minCV) # find good regions for the Good read
        sys.stderr.write("covered bps: {}\naverage coverage depth: {}\n".format(cov_bps, avg_cov_depth))
        if len(good_regions) == 0 :
            sys.stderr.write("PacBio read does not have any good region covered by the Illumina reads")
        else: # examine good regions one by one
            # first print out all good regions
            sys.stderr.write("Good regions:\n")
            for i in good_regions:
                sys.stderr.write('({}, {}): {}\t'.format(i[0], i[1], i[1] - i[0])) # (start_pos, end_pos): length of this good region
            sys.stderr.write('\n')

            seqName = os.path.basename(args.seqFile).split('.')[0] # e.g. m130828_041445_00123_c100564312550000001823090912221381_s1_p0__58103__7045_8127.fasta
            seqName = re.sub('__','/',seqName) # change __ back to /
            # try to correct PacBio sequence at each good region
            for good_region_index in xrange(len(good_regions)):
                sys.stderr.write("====\nworking on region {}\n".format(good_region_index))
                # step 1 - find consensus, polymorphic positions, and coverage depths for the PacBio read
                poly_bps, poly_ins, consensus_bps, consensus_ins, cvs = metalrec_lib.get_poly_pos(ref_bps, ref_ins_dict, good_regions[good_region_index])
                # step 2 - extend the PacBio sequence to include the insertion positions, and find the correspondance between positions from original and extened sequences
                newSeq, bp_pos_dict, ins_pos_dict = metalrec_lib.ref_extension(poly_bps, poly_ins, consensus_bps, consensus_ins, rseq, region=good_regions[good_region_index],print_width=args.width, verbose=args.verbose)
                # step 3 - update consensus and polymorphic positions according to the new positons in the extended sequence
                poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext = metalrec_lib.update_pos_info(poly_bps, poly_ins, consensus_bps, consensus_ins, bp_pos_dict, ins_pos_dict)
                # step 4 - make array to indicate the type for each position on the extended PacBio sequence
                type_array, coordinates = metalrec_lib.make_type_array(poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext,verbose=args.verbose)
                # step 5 - construct array for all the reads that passed the specified threshold, and number of repeats for each unique read (single or paired)
                read_array, read_counts = metalrec_lib.make_read_array(read_info, bp_pos_dict, ins_pos_dict, type_array, poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext, ext_region=coordinates)
                if not numpy.any(read_array != 0): # if all read calls for 0
                    sys.stderr.write("PacBio read does not have any good reads covering the region.\n")
                else:
                    refOut = open(args.oSeqFile, 'a') # output file for the corrected PacBio sequence (contigs if it is split)
                    #shortOut = open(shortSeqFile,'a')

                    # step 6 - find error corrected sequence by filling the gaps in the greedy fashion
                    if args.verbose:
                        region_outdir = args.outDir+'region' + str(good_region_index)
                        fastaFile = args.outDir + 'goodreads.fasta'
                    else:
                        region_outdir = None
                        fastaFile = None
                    ref_new = metalrec_lib.fill_gap(read_array,args.minOverlap,args.minOverlapRatio, fastaFile, region_outdir, read_info, verbose=args.verbose)
                    # step 7 - convert the array for the new PacBio sequence to string of nucleotides
                    contiguous_seqs = metalrec_lib.split_at_gap(ref_new[2], ref_new[0])
                    contiguous_lengths = numpy.array(map(len, contiguous_seqs))
                    max_ind = numpy.argmax(contiguous_lengths)
                    # in verbose mode, print the comparison between the original sequence, the extended sequence, and the corrected sequence
                    ## write the newly corrected sequence to the output sequence file
                    # header format: >1 (0, 1048) gap length: 16
                    #refOut.write('>{}/{}_{}_M ({}, {}) length: {}\n{}\n'.format(seqName, good_region_index,max_ind, good_regions[good_region_index][0], good_regions[good_region_index][1], ref_new[1], contiguous_seqs[max_ind]))
                    refOut.write('>{}/{} ({}, {}) scrub; length: {}\n'.format(seqName, good_region_index, good_regions[good_region_index][0], good_regions[good_region_index][1], ref_new[1]))
                    start = 0
                    outString = contiguous_seqs[max_ind]
                    while start < contiguous_lengths[max_ind]:
                        refOut.write(outString[start:min(start+100, len(outString))] + "\n")
                        start = start + 100
                    #if len(contiguous_seqs) > 1:
                    #    for i in xrange(len(contiguous_seqs)):
                    #        if i != max_ind:
                    #            shortOut.write('>{}/{}_{} ({}, {}) length: {}\n{}\n'.format(seqName, good_region_index, i,  good_regions[good_region_index][0], good_regions[good_region_index][1], contiguous_lengths[i], contiguous_seqs[i]))

                    refOut.close()
                    #shortOut.close()
    sys.stderr.write("total time :" + str(time.time() - start_time) +  " seconds")
    sys.stderr.write("\n===========================================================\nDone\n")
##==============================================================
## call from command line (instead of interactively)
##==============================================================

if __name__ == '__main__':
    sys.exit(main())
