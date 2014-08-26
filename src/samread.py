#!/usr/bin/python

import metalrec_lib
import re
import sys
revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1]) # find reverse complement of a DNA sequence
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
        self.mate = None

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
    def is_record_bad(self, maxSub=3, maxIns=3, maxDel=3,maxSubRate=0.02, maxInsRate=0.2, maxDelRate=0.2):
        return metalrec_lib.is_record_bad(self.alignRecord, maxSub, maxIns, maxDel,maxSubRate, maxInsRate, maxDelRate)

    # get nucleotide, base by base
    def get_bases(self):
        cigar_string = self.cigarstring
        qseq = self.qSeq
        start_pos = self.rstart
        return metalrec_lib.get_bases(self.cigarstring, self.qSeq, self.rstart)
    
    # trim the query sequence to only keep the aligned part, get rid of the clipped sequence at the two ends
    # only soft clipped part is present in the qseq, not the hard clipped sequence
    def get_trim_qseq(self):
        if 'S' in self.cigarstring:
            cigar_list = re.findall('\D|\d+',self.cigarstring) # cigar string as a list of str (numbers and operations) 
            if cigar_list[1] == 'S':
                qSeq = self.qSeq[int(cigar_list[0]):]
            if cigar_list[-1] == 'S':
                qSeq = self.qSeq[:int(cigar_list[-2])]
            return qSeq
        else:
            return self.qSeq

    # get the trimmed original read sequence
    def get_read_seq(self):
        read_seq = get_trim_qseq(self)
        if is_reverse(self):
            return revcompl(read_seq)
        else:
            return read_seq

    # generate the alignment record for the read, also considering its mate in this function
    # record will be written if a read is mapped well, or at least one of the read of a pair is well mapped
    def generate_sam_record(self):
        # [0] qname stays the same

        # [1] new flag: mapped/unmapped status might change
        if not is_unmapped(self) and is_record_bad(self, maxSub=3, maxIns=3, maxDel=3,maxSubRate=0.02, maxInsRate=0.2, maxDelRate=0.2): # read was originally mapped but didn't pass the threshold
            self.flag += 0x4 # change its own "unmapped" flag
            if self.mate is not None:
                self.mate.flag += 0x8 # change mate's "mate_unmapped" flag
                self.cigarstring = '*'
        if not mate_is_unmapped(self) and is_record_bad(self.mate, maxSub=3, maxIns=3, maxDel=3,maxSubRate=0.02, maxInsRate=0.2, maxDelRate=0.2): # mate was originally mapped but didn't pass the threshold
            self.mate.flag += 0x8 # if mate record is bad, mark its mate as unmapped
            if self.mate is not None:
                self.flag += 0x4 # if mate record is bad, mark its mate as unmapped
                self.mate.cigarstring = '*'
        self.fields[1] = str(self.flag)
        if self.mate is not None:
            self.mate.fields[1] = str(self.mate.flag)

        # [2] rname stays the same

        # [3] pos: 1-based leftmost mapping POSition might change after re-mapping 
        self.fields[3] = self.rstart
        if self.mate is not None:
            self.mate.fields[3] = self.mate.rstart

        # [4] MAPQ stays the same

        # [5] CIGAR string might have changed after re-mapping to PacBio sequence
        self.fields[5] = self.cigarstring
        if self.mate is not None:
            self.mate.fields[5] = self.mate.cigarstring

        # [6] RNEXT Ref. name of the mate/next read: should always be = since there is only 1 PacBio sequence 

        # [7] PNEXT Position of the mate/next read: might change because of re-mapping
        self.fields[7] = self.mate.rstart
        if self.mate is not None:
            self.mate.fields[7] = self.rstart

        # [8] SEQ segment SEQuence: trimmed original sequence (only the part mapped to PacBio sequence)
        self.fields[8] = get_read_seq(self)
        if self.mate is not None:
            self.mate.fields[8] = get_read_seq(self.mate)

        # [9] QUAL stays the same

        ## Paste all fields together, separated by tabs, if at least one read of the pair is mapped. 
        ## Two lines if read has a mate, no matter if they are mapped or not
        if not is_unmapped(self) or (self.mate is not None and not mate_is_unmapped(self)):
            record = '\t'.join(self.fields)
            if self.mate is not None:
                record += '\n' + '\t'.join(self.mate.fields)
        else:
            record = ''
        return record
