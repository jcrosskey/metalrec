#!/usr/bin/python

import metalrec_lib
import re
import sys
''' Class SamRead '''

class SamRead:
    ''' Class representing an aligned read, particularly from sam file, instead of bam file '''
    def __init__(self, alignRecord):
        self.alignRecord = alignRecord.strip('\n')
        fields = alignRecord.split('\t')
        self.fields = fields
        self.qname = fields[0] # Query template name
        self.flag = int(fields[1]) # bitwise Flag
        self.rname = fields[2] # reference sequence name
        self.rstart = int(fields[3]) # starting mapping position on the reference sequence, 1-based
        try:
            self.mapQ = int(fields[4]) # mapping quality
        except ValueError: # if map value is not integer
            self.mapQ = -1
        self.cigarstring = fields[5] # CIGAR string
        self.mate_rname = fields[6] # reference for mate/next read
        self.mate_rstart = fields[7] # position of the mate/next read
        self.tLen = abs(int(fields[8]))  # observed template length. If all segments are mapped to the same reference, 
                                         # the unsinged observed template length equals the number of bases drom the leftmost mapped base to the right most mapped base.
                                         # The leftmost segment has a plus sign and the rightmost has a minus sign.
                                         # In other words, it's the insert size
        self.qSeq = fields[9] # segment sequence

    def get_tags(self):
        tags = dict() # other tags available
        for i in xrange(10,len(self.fields)):
            tag_name = self.fields[i][:2]
            tag_value = self.fields[i][5:]
            tags[tag_name] = tag_value
        return tags

    def get_cigar_info(self):
        return metalrec_lib.cigar(self.cigarstring)

    def get_rend(self):
        return self.rstart + self.get_cigar_info()['ref_len'] - 1 # ending mapping position on the reference sequence, 1-based

    # FLAG explanation one by one

    def is_paired(self):
        return self.flag & 0x1 == 0x1 # 0x1: read is paired

    def is_proper_pair(self):
        return self.flag & 0x2 == 0x2 # 0x2: reads mapped in proper pair 

    def is_unmapped(self):
        return self.flag & 0x4 == 0x4 # 0x4: read unmapped

    def mate_is_unmapped(self):
        return self.flag & 0x8 == 0x8 # 0x8: mate/next read unmapped

    def is_reverse(self):
        return self.flag & 0x10 == 0x10 # 0x10: read being reverse complemented

    def mate_is_reverse(self):
        return self.flag & 0x20 == 0x20 # 0x20: mate/next read being reverse complemented

    def is_read1(self):
        return self.flag & 0x40 == 0x40 # 0x40: first read in the pair

    def is_read2(self):
        return self.flag & 0x80 == 0x80 # 0x80: second read in the pair

    def is_secondary(self):
        return self.flag & 0x100 == 0x100 # 0x100: secondary alignment, false if mapping is primary

    def is_qcfail(self):
        return self.flag & 0x200 == 0x200 # 0x200: not passing quality control

    def is_duplicate(self):
        return self.flag & 0x400 == 0x400 # 0x400: read is PCR or optical duplicate

    def is_supplementary(self):
        return self.flag & 0x800 == 0x800 # 0x800: supplementary alignment (part of a chimeric alignment)
    
    # check and see if this mapping is too noisy to be included
    def is_record_bad(self, maxSub=3, maxIns=3, maxDel=3,maxErrRate=0.20):
        return metalrec_lib.is_record_bad(self.alignRecord, maxSub, maxIns, maxDel, maxErrRate)

    # get nucleotide, base by base
    def get_bases(self):
        ''' from CIGAR string, query segment (in alignment record), and starting position (1-based) on the ref sequence, return position wise base call from the read.
            Assuming that cigar_string is available, doesn't matter if it's sam 1.3 format or 1.4 format
            
            Input: cigar_string, query aligned segment, and the 1-based starting mapping position on the reference sequence
            
            Output: (pos_dict, ins_dict), tuple of 2 dictionaries, one for the non-insertion positions, and one for the insertion positions.
                    When there is a deletion from the reference sequence, the base called will be "D"
        '''
        cigar_string = self.cigarstring
        qseq = self.qSeq
        start_pos = self.rstart
        return metalrec_lib.get_bases(self.cigarstring, self.qSeq, self.rstart)
