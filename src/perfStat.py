#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
perfStat

Created on Fri May 23 13:14:18 2014

@author: cjg
"""
import sys
import re
import argparse
## =================================================================
## trim_name
## =================================================================
def trim_name(name):
    '''
    trim a read's name so that it corresponds to the name of the filtered subread from which the sequence was generated.
    '''
    m = re.search('\/\d+_\d+',name)
    if m is None:
        sys.stdout.write("read name {} does not have the right format\n".format(name))
        return name
    else:
        return m.string[:m.end()]
## =================================================================
## perfstat
## =================================================================
def perfStat(statFile, outFile, delta=0):
    '''
    Find the statistics for the perfectly corrected (mapped to reference genomes) ZMWs.
    Input:  statFile - stat.read (output from samStat.py)
    Output: outFile - as in the description
    '''
    perf_dict = dict()
    with open(statFile, 'r') as inStat:
        inStat.readline()
        for line in inStat:
            line = line.strip("\n").split("\t")
            readName = line[0]
            mapRatio = float(line[5])
            readLen = int(line[3])
            nEdit = int(line[12])
            Mapping = int(line[1])
            if Mapping == 0 and abs(mapRatio-1) <= delta: # perfect primary mapping
                FR = trim_name(readName) # filtered subread name
                if FR not in perf_dict:
                    perf_dict[FR] = {'region':0, 'maxLen':0, 'totLen':0, 'mapRatio':0}
                perf_dict[FR]['region'] += 1
                if perf_dict[FR]['maxLen'] < readLen:
                    perf_dict[FR]['maxLen'] =  readLen
                    perf_dict[FR]['mapRatio'] = mapRatio
                perf_dict[FR]['totLen'] += readLen
    outStat = open(outFile, 'w')
    outStat.write('readName\tregions\tmaxLen\tmaxLenRatio\ttotLen\n')
    for key in sorted(perf_dict):
        d = perf_dict[key]
        outStat.write('{}\t{}\t{}\t{}\t{}\n'.format(key, d['region'], d['maxLen'],d['mapRatio'], d['totLen']))
    outStat.close()
## =================================================================
## argument parser
## =================================================================
parser = argparse.ArgumentParser(description="Find the statistics for the perfectly corrected PB reads",
                                 prog = 'samStat', #program name
                                 prefix_chars='-', # prefix for options
                                 fromfile_prefix_chars='@', # if options are read from file, '@args.txt'
                                 conflict_handler='resolve', # for handling conflict options
                                 add_help=True, # include help in the options
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter # print default values for options in help message
                                 )
## input files and directories
parser.add_argument("-i","--in",help="input read mapping stat file",dest='statFile',required=True)
## output directory
parser.add_argument("-o","--out",help="output statistics summary file",dest='outputFile',required=True)
## option
parser.add_argument("-d","--delta",help="tolerance of mapping quality, deviation from 1",dest='delta',type=float, default=0.0)


## =================================================================
## main function
## =================================================================
def main(argv=None):
    if argv is None:
        args = parser.parse_args()

    perfStat(args.statFile, args.outputFile, args.delta)

##==============================================================
## call from command line (instead of interactively)
##==============================================================

if __name__ == '__main__':
    sys.exit(main())
