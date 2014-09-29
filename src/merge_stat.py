#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Merge the stat.read files from 2 different mappings and compare the result
Created on Tue Sep  9 11:32:48 EDT 2014
@author: cjg
"""
import sys
import os
import argparse
import time
import re

## =================================================================
## argument parser
## =================================================================
parser = argparse.ArgumentParser(description=" Merge the stat.read files from 2 different mappings and compare the result",
                                 prog = 'merge_stat', #program name
                                 prefix_chars='-', # prefix for options
                                 fromfile_prefix_chars='@', # if options are read from file, '@args.txt'
                                 conflict_handler='resolve', # for handling conflict options
                                 add_help=True, # include help in the options
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter # print default values for options in help message
                                 )

## input files and directories
parser.add_argument("statFiles",help="2 stat files to merge: before and after Error Correction",nargs=2)
## output directory
parser.add_argument("output",help="output file", nargs='?', type=argparse.FileType('w'), default=sys.stdout)

def read_EC_stat(read_stat_file):
    read_dict = dict()
    with open(read_stat_file,'r') as stat:
        header = stat.readline()
        for line in stat:
            readName = line.strip("\n").split("\t")[0]
            if readName.count(':') == 2:
                sName = readName.split(":")[0]
            else:
                sName = '/'.join(readName.split("/")[:-1])
            print sName
            readLen = map(int,re.findall('\d+',readName[readName.find('(') : readName.find(')')]))
            readLen = readLen[1] - readLen[0]
            nMatchBp, nInsBp, nDelBp, nSubBp, nEdit = map(int, line.strip("\n").split("\t")[2:7])
            if sName not in read_dict:
                read_dict[sName] = [readLen, nMatchBp, nInsBp, nDelBp, nSubBp, nEdit, 1]
            else:
                read_dict[sName] = [read_dict[sName][i] + [readLen, nMatchBp, nInsBp, nDelBp, nSubBp, nEdit, 1][i] for i in xrange(7)]
    return read_dict

def extract_from_allFR(read_stat_file, read_dict, out):
    out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format("readName", "EC_readLen", "EC_nMatchBp", "EC_nInsBp", "EC_nDelBp", "EC_nSubBp", "EC_nEdit", "EC_regions", "readLen", "nMatchBp", "nInsBp", "nDelBp", "nSubBp", "nEdit"))
    with open(read_stat_file, 'r') as stat:
        header = stat.readline()
        for line in stat:
            line = line.strip("\n").split("\t")
            readName = line[0]
            readLen = readName.split("/")[-1]
            readLen = map(int, readLen.split("_"))
            readLen = readLen[1] - readLen[0]
            if readName in read_dict:
                read_dict[readName] += [readLen, line[2], line[3], line[4], line[5], line[6]]
    for read in sorted(read_dict):
        out.write('{}\t{}\n'.format(read, "\t".join(map(str,read_dict[read]))))



def get_length(lengthFile, trim=""):
    length_dict = dict()
    with open(lengthFile, 'r') as length:
        for line in length:
            seqName, seqLen = line.strip('\n').split("\t")
            if trim != '':
                seqName = seqName.split(trim)[0]
            length_dict[seqName] = int(seqLen)
    return length_dict

def read_stat(statFile, splitChar='',length_dict=dict()):
    ''' scan stat.read file and store the information in dictionary
    '''
    mydict = dict()
    with open(statFile, 'r') as stat:
        stat.readline()
        for line in stat:
            line = line.strip("\n")
            readName, nMapping, Mappings = line.split("\t")
            if splitChar != '':
                readName_base = splitChar.join(readName.split(splitChar)[:-1]) # if there is split char, 
                coord = readName.split(splitChar)[-1]
                coord = re.findall('\d+', coord)
                if len(coord) == 1: # only :0 :1 , etc
                    if readName in length_dict:
                        segLen = length_dict[readName]
                    else:
                        segLen = 'NA'
                elif len(coord) == 2: # e.g. 0_1284
                    coord = map(int, coord)
                    segLen = coord[1] - coord[0]
            else:
                segLen = 'NA'
            mappings = Mappings.split(")(")
            for mapping in mappings:
                if mapping.split('#')[2].strip() == '1': # primary alignment
                    at_pos = mapping.index('@')
                    at_pos_end = mapping.index(',',at_pos)
                    map_positions = map(int,re.findall('\d+', mapping.split('#')[0]))
                    mapLen = map_positions[1] - map_positions[0] + 1
                    NM = int(mapping[ (at_pos + 2) : at_pos_end ] )
                    break
            if readName_base not in mydict:
                mydict[readName_base] = [[nMapping, segLen, mapLen, NM, Mappings]]
            else:
                mydict[readName_base].append([nMapping, segLen, mapLen, NM, Mappings])
    return mydict


## =================================================================
## main function
## =================================================================
def main(argv=None):
    if argv is None:
        args = parser.parse_args()
    statFile1, statFile2 = args.statFiles # input stat files
    myout = args.output # output file object
    read_dict = read_EC_stat(statFile2)
    extract_from_allFR(statFile1, read_dict, myout)

##==============================================================
## call from command line (instead of interactively)
##==============================================================

if __name__ == '__main__':
    sys.exit(main())

