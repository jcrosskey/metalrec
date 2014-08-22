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
parser.add_argument("-o","--out",help="output corrected PacBio sequence file",dest='oSeqFile',required=True)
parser.add_argument("-rs","--redSam",help="reduced sam file with only good alignment",dest='redSam',required=False, default = '') # optional, for visualization
#parser.add_argument("-os","--outSam",help="alignment of same Illumina sequences to the corrected PacBio read",dest='outSam',required=True)
parser.add_argument("-of","--outFasta",help="fasta file including reads that are reserved from the mapping results",dest='outFasta',required=True)
parser.add_argument("-od","--outDir",help="directory for the intermediate files",dest='outDir',default = None)

# options
parser.add_argument("-m",help="read mode, s(ingle) or p(air)",dest='rMode',default='p',choices=['s','p'])
parser.add_argument("-v","--verbose",help="verbose, more output",action='store_true',dest='verbose')

## setting thresholds
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
## =================================================================
## main function
## =================================================================
def main(argv=None):
    # parse command line arguments    
    if argv is None:
        args = parser.parse_args()

    # check input and output file settings
    if not os.path.exists(args.seqFile):
        sys.exit("input PacBio sequence file does not exist!\n")
    if not os.path.exists(args.samFile):
        sys.exit("input sam file does not exist!\n")
    if os.path.exists(args.oSeqFile):
        os.remove(args.oSeqFile)
        sys.stdout.write("Output sequence file already exists, overwrite.\n")

    # read the PacBio sequence into memory
    rseq = metalrec_lib.read_single_seq(args.seqFile)
    # process sam file and save the read info, depending on the reading mode 
    if args.rMode == 'p':
        ref_bps, ref_ins_dict, read_info = metalrec_lib.read_and_process_sam_pair(args.samFile, rseq, maxSub=args.maxSub, maxIns=args.maxIns, maxDel=args.maxDel,maxSubRate=args.maxSubRate, maxInsRate=args.maxInsRate, maxDelRate=args.maxDelRate, minPacBioLen=args.minPacBioLen, minCV=args.minCV, outsam=args.redSam, outFastaFile=args.outFasta)
    else:
        ref_bps, ref_ins_dict, read_info = metalrec_lib.read_and_process_sam(args.samFile, rseq, maxSub=args.maxSub, maxIns=args.maxIns, maxDel=args.maxDel,maxSubRate=args.maxSubRate, maxInsRate=args.maxInsRate, maxDelRate=args.maxDelRate, minPacBioLen=args.minPacBioLen, minCV=args.minCV, outsam=args.redSam, outFastaFile=args.outFasta)

    good_regions = metalrec_lib.get_good_regions(ref_bps, rseq, minPacBioLen=args.minPacBioLen, minCV=args.minCV) # find good regions for the PacBio read
    print good_regions
    if len(good_regions) == 0 :
        sys.exit("PacBio read does not have good region covered by the Illumina reads")
    else: # examine good regions one by one
        refOut = open(args.oSeqFile, 'a')
        for good_region_index in xrange(len(good_regions)):
            # step 1 - find consensus, polymorphic positions, and coverage depths for the PacBio read
            poly_bps, poly_ins, consensus_bps, consensus_ins, cvs = metalrec_lib.get_poly_pos(ref_bps, ref_ins_dict, good_regions[good_region_index])
            # step 2 - extend the PacBio sequence to include the insertion positions, and find the correspondance between positions from original and extened sequences
            newSeq, bp_pos_dict, ins_pos_dict = metalrec_lib.ref_extension(poly_bps, poly_ins, consensus_bps, consensus_ins, rseq)
            # step 3 - update consensus and polymorphic positions according to the new positons in the extended sequence
            poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext = metalrec_lib.update_pos_info(poly_bps, poly_ins, consensus_bps, consensus_ins, bp_pos_dict, ins_pos_dict)
            # step 4 - make array to indicate the type for each position on the extended PacBio sequence
            type_array = metalrec_lib.make_type_array(poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext)
            # step 5 - construct array for all the reads that passed the specified threshold, and number of repeats for each unique read (single or paired)
            read_array, read_counts = metalrec_lib.make_read_array(read_info, bp_pos_dict, ins_pos_dict, type_array, poly_bps, poly_ins, consensus_bps, consensus_ins)
            # step 6 - find error corrected sequence by filling the gaps in the greedy fashion
            ref_new, gap = metalrec_lib.fill_gap(read_array, args.outFasta, args.outDir, read_info)
            # step 7 - convert the array for the new PacBio sequence to string of nucleotides
            ref_new_short, ref_new_long = metalrec_lib.array_to_seq(ref_new)
            # in verbose mode, print the comparison between the original sequence, the extended sequence, and the corrected sequence
            if args.verbose:
                pass # need to fill in this part later TODO
            ## write the newly corrected sequence to the output sequence file
            # header format: >1 (0, 1048) gap length: 16
            refOut.write('>{} ({}, {}) gap length: {}\n{}\n'.format(good_region_index, good_regions[good_region_index][0], good_regions[good_region_index][1], gap, ref_new_short))
        refOut.close()
##==============================================================
## call from command line (instead of interactively)
##==============================================================

if __name__ == '__main__':
    sys.exit(main())
