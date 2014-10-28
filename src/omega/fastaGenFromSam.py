#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
samStat

Created on Fri May 23 13:14:18 2014

@author: cjg
"""
import sys
import argparse
import pysam
import samread

## =================================================================
## samStat function with SamRead class
## =================================================================
def samStat(samFile, outputFile):
    ''' From resulted sam or bam file of mapping, find information of reference sequences and reads.
        For reference sequences: 
        1. coverage percentage
        2. coverage depth at each base pair
        3. error rate/number (ins, del, sub)
        4. number of reads mapped to it
        
        For reads:
        1. number of reported alignments that contains the query read
        2. for each such alignment, what's the reference name and the qregion of the read
        
        Input:
        1. samFile: sam (bam) file name
        2. outputFile: file for writing output
        3. fileformat: sam for now (should be either sam or bam, should do auto detect..)
    '''
    nReferences = 0 # number of reference sequences
    refLens = [] # list of reference length
    refNames = []# list of reference names
    count = 0 # number of aligned records in the sam file, could be either mapped or unmapped

    # dictionaries for the reference sequences and the read sequences
    refSeq_dict = dict()
    readSeq_dict = dict()

    sys.stdout.write(">> Scan sam file \n")
    # start scanning sam file
    with open(samFile,'r') as mysam:
        for line in mysam:
            if line[0] == '@': # header line
                if line[1:3] == 'SQ': # reference sequence dictionary
                    nReferences += 1
                    rname = line[(line.find('SN:') + len('SN:')) : line.find('\t',line.find('SN:'))] # referenece sequence name
                    rLen = line[(line.find('LN:') + len('LN:')) : line.find('\t',line.find('LN:'))] # reference sequence length
                    refLens.append(int(rLen))
                    refNames.append(rname)
            else: # non-header line
                line = line.strip()            
                count += 1 # number of read record increased by 1
                myread = samread.SamRead(line) # parse the alignment record, using SamRead class 

                if myread.cigarstring == '*' or myread.is_unmapped(): # read is not mapped
                    continue
                
                # if this reference sequence is not in the dictionary, initiate it
                if not refSeq_dict.has_key(myread.rname):
                    refLen = refLens[refNames.index(myread.rname)] # length of the reference sequence
                    refSeq_dict[myread.rname] = {'refLen':refLen, 'nReads':0, 'nReadsBp':0, 'nMatchBp':0,'nInsBp':0, 'nDelBp':0, 'nSubBp':0, 'nEdit':0,'coverage':[0]*refLen}

                if not readSeq_dict.has_key(myread.qname):
                    readSeq_dict[myread.qname] = {'nMapping':0, 'mapInfo':[]}

                #print qname, '\t', rname, '\t', refLen

                ## update the dictionary corresponding to the reference sequence
                refSeq_dict[myread.rname]['nReads'] += 1 # update number of mapped reads
                readSeq_dict[myread.qname]['nMapping'] += 1 # update number of mappings

                ## get CIGAR info
                cigarInfo = myread.get_cigar_info()

                # update number of bps mapped to this ref seq
                refSeq_dict[myread.rname]['nReadsBp'] += cigarInfo['ref_len'] 
                # update matching and substitution bps if possible
                if cigarInfo['match_len'] is not None:
                    refSeq_dict[myread.rname]['nMatchBp'] += cigarInfo['match_len']
                if cigarInfo['sub_len'] is not None:
                    refSeq_dict[myread.rname]['nSubBp'] += cigarInfo['sub_len']

                refSeq_dict[myread.rname]['nInsBp'] += cigarInfo['ins_len'] # update number of insertion bps
                refSeq_dict[myread.rname]['nDelBp'] += cigarInfo['del_len'] # update number of deletion bps

                # update edit distance
                if 'NM' in myread.get_tags():
                    myread.NM = int(myread.get_tags()['NM'])
                    refSeq_dict[myread.rname]['nEdit'] += myread.NM

                # update the coverage at the mapped positions
                for apos in xrange(myread.rstart-1, myread.get_rend()):
                    refSeq_dict[myread.rname]['coverage'][apos] += 1

                # store the mapping information for this read:
                # start and end positions for both the query read and the ref seq
                # is this a secondary alignment?
                # is this a reverse complement?
                readSeq_dict[myread.qname]['mapInfo'].append([cigarInfo['match_len'], cigarInfo['ins_len'], cigarInfo['del_len'], cigarInfo['sub_len'], myread.NM, myread.rstart, myread.get_rend(), myread.is_secondary(), myread.is_reverse(),myread.rname])

                if count % 10000 == 0:
                    sys.stdout.write('  scanned {} records\n'.format(count))
                
    sys.stdout.write(">> Write statistics in output file \n")

    # print out statistics information for the reference sequences
    myout1 = open(outputFile+".ref", 'w')
    myout1.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format('refName', 'refLen','nReads', 'nReadsBp', 'nMatchBp','nInsBp','nDelBp','nSubBp','nEdit','nCovBp','maxCov','avgCov','coverage'))

    for key in refSeq_dict:
        d = refSeq_dict[key]
        nCovBp = d['refLen'] - d['coverage'].count(0) # number of covered base pairs = total bps - 0-coverage bps
        myout1.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(key, d['refLen'],d['nReads'], d['nReadsBp'], d['nMatchBp'],d['nInsBp'],d['nDelBp'],d['nSubBp'],d['nEdit'],nCovBp,max(d['coverage']),float(d['nReadsBp'])/float(d['refLen']),float(nCovBp)/float(d['refLen'])))

    myout1.close()

    # print out statistics information for the reads
    myout2 = open(outputFile+".read", 'w')
    myout2.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format('readName', 'Mappings','nMatchBp', 'nInsBp', 'nDelBp', 'nSubBp', 'nEdit', 'rstart','rend', 'is_secondary','is_revcomp','rname'))
    for key in sorted(readSeq_dict):
        d = readSeq_dict[key]
        #myout2.write("{}\t{}\t".format(key, d['nMapping']))
        for i in xrange(d['nMapping']):
            thismap = d['mapInfo'][i]
            # qstart, qend # rstart, rend # secondary # forward/backward @  edit distance, refName
            myout2.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(key,i,thismap[0],thismap[1],thismap[2],thismap[3],thismap[4], thismap[5], thismap[6], thismap[7], thismap[8], thismap[9].split()[0]))

    myout2.close()
    return readSeq_dict
## =================================================================
## argument parser
## =================================================================
parser = argparse.ArgumentParser(description="parse sam file and get summary statistics",
                                 prog = 'samStat', #program name
                                 prefix_chars='-', # prefix for options
                                 fromfile_prefix_chars='@', # if options are read from file, '@args.txt'
                                 conflict_handler='resolve', # for handling conflict options
                                 add_help=True, # include help in the options
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter # print default values for options in help message
                                 )

## input files and directories
parser.add_argument("-i","--in",help="input sam file",dest='samFile',required=True)

## output directory
parser.add_argument("-o","--out",help="output statistics file",dest='outputFile',required=True)

## =================================================================
## main function
## =================================================================
def main(argv=None):
    
    if argv is None:
        args = parser.parse_args()

    readSeq_dict = samStat(args.samFile,args.outputFile)

##==============================================================
## call from command line (instead of interactively)
##==============================================================

if __name__ == '__main__':
    sys.exit(main())
