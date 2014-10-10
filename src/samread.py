#!/usr/bin/python

import metalrec_lib
import re
import sys
sys.path.append("/lustre/atlas/scratch/chaij1/csc124/biopython-1.64/lib/python2.7/site-packages")
import math
from Bio import pairwise2 # pairwise alignment using dynamic programming
from Bio.pairwise2 import format_alignment

''' reverse complement function '''
revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}.get(B,'N') for B in x ][::-1]) # find reverse complement of a DNA sequence


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
    
    # check and see if this mapping is too noisy to be included, if so, change the flag indicating if the read was mapped
    def is_read_bad(self, maxSub=-1, maxIns=-1, maxDel=-1,maxSubRate=0.05, maxInDelRate=0.3):
        is_bad = metalrec_lib.is_record_bad(self.alignRecord, maxSub, maxIns, maxDel,maxSubRate, maxInDelRate)
        if is_bad and not self.is_unmapped():
            self.flag += 0x4 # change its own "unmapped" flag
        return is_bad

    # get nucleotide, base by base
    def get_bases(self):
        cigar_string = self.cigarstring
        qseq = self.qSeq
        start_pos = self.rstart
        return metalrec_lib.get_bases(self.cigarstring, self.qSeq, self.rstart)
    
    # trim the query sequence to only keep the aligned part, get rid of the clipped sequence at the two ends
    # only soft clipped part is present in the qseq, not the hard clipped sequence
    def trim_qseq(self):
        if 'S' in self.cigarstring:
            cigar_list = re.findall('\D|\d+',self.cigarstring) # cigar string as a list of str (numbers and operations) 
            if cigar_list[1] == 'S':
                self.qSeq = self.qSeq[int(cigar_list[0]):]
                del cigar_list[:2]
            if cigar_list[-1] == 'S':
                self.qSeq = self.qSeq[:-int(cigar_list[-2])]
                del cigar_list[-2:]
            self.cigarstring = ''.join(cigar_list)

    ## Retrieve the clipped part as long as it's still inside the long PacBio sequence.
    def get_seg(self, rseq):
        ''' This function does not alter the qSeq of the read. Get the segment for realigning to PacBio sequence, considering the clipped part.
            Input:  rseq - reference sequence (PacBio)
            Output: keep_left, keep_right - number of clipped bases kept to the left and right of the mapped segment of Illumina read
                    seg - segment that totally reside inside the PacBio read
        '''
        keep_left = 0
        keep_right = 0
        seg = str(self.qSeq)
        if 'S' in self.cigarstring: # if it's hard clipping, the sequence information is not available in sam file
            cigar_info = metalrec_lib.cigar(self.cigarstring)
            #print "before clipping: ", seg #DEBUG
            #print "starting position on the ref sequence is ", self.rstart
            #print "left clip: ", cigar_info['left_clip_len']
            #print "right clip: ", cigar_info['right_clip_len']
            if cigar_info['left_clip_len'] > 0 or cigar_info['right_clip_len'] > 0:
                if cigar_info['left_clip_len'] > 0:
                    keep_left = min( self.rstart - 1, cigar_info['left_clip_len'])
                    #print "keep_left: ", keep_left
                    seg = seg[(cigar_info['left_clip_len'] - keep_left):] # change the segment sequence, instead of the whole read sequence
                if cigar_info['right_clip_len'] > 0:
                    keep_right = min(len(rseq) - self.get_rend(), cigar_info['right_clip_len'])
                    #print "keep_right: ", keep_right
                    if keep_right < cigar_info['right_clip_len']:
                        seg = seg[:-(cigar_info['right_clip_len'] - keep_right)]
        #print "after clipping: ", seg
        return keep_left, keep_right, seg
    ## trim query sequences, with clipped part in consideration. This function does not change cigar string.
    ## It takes care of soft clipping, not hard clipping
    ## Retrieve the clipped part as long as it's still inside the long PacBio sequence.
    def trim_qseq_clip(self, rseq):
        ''' This function alters the qSeq of the read. Get the segment for realigning to PacBio sequence, considering the clipped part.
            Input:  rseq - reference sequence (PacBio)
            Output: keep_left, keep_right - number of clipped bases kept to the left and right of the mapped segment of Illumina read
        '''
        keep_left = 0
        keep_right = 0
        if 'S' in self.cigarstring: # if it's hard clipping, the sequence information is not available in sam file
            cigar_info = metalrec_lib.cigar(self.cigarstring)
            #print "before clipping: ", self.seg #DEBUG
            #print "starting position on the ref sequence is ", self.rstart
            #print "left clip: ", cigar_info['left_clip_len']
            #print "right clip: ", cigar_info['right_clip_len']
            if cigar_info['left_clip_len'] > 0 or cigar_info['right_clip_len'] > 0:
                if cigar_info['left_clip_len'] > 0:
                    keep_left = min( self.rstart - 1, cigar_info['left_clip_len'])
                    #print "keep_left: ", keep_left
                    self.seg = self.seg[(cigar_info['left_clip_len'] - keep_left):] # change the segment sequence, instead of the whole read sequence
                if cigar_info['right_clip_len'] > 0:
                    keep_right = min(len(rseq) - self.get_rend(), cigar_info['right_clip_len'])
                    #print "keep_right: ", keep_right
                    if keep_right < cigar_info['right_clip_len']:
                        self.seg = self.seg[:-(cigar_info['right_clip_len'] - keep_right)]
        #print "after clipping: ", self.seg
        return keep_left, keep_right

    # get the read segment, could be hard clipped or not, depending on if the clip function has been called
    def get_read_seq(self, keep_orientation=False):
        ''' Get the read segment.
            Input:  keep_orientation - switch to keep the original orientation of the read, or set it the same as how it's aligned
                                       default: keep orientation as how it's mapped to the PacBio sequence
            Output: read_seq - read segment (or the whole read)
        '''
        read_seq = self.qSeq
        if self.is_reverse() and keep_orientation:
            return revcompl(read_seq)
        else:
            return read_seq

    # generate the alignment record for the read, also considering its mate in this function
    # record will be written if a read is mapped well, or at least one of the read of a pair is well mapped
    def generate_sam_record(self, maxSub=-1, maxIns=-1, maxDel=-1,maxSubRate=0.1, maxInDelRate=0.3):
        # [0] qname stays the same
        # [1] new flag: mapped/unmapped status might change
        if (not self.is_unmapped()) and (self.mate is not None) and (self.is_read_bad(maxSub, maxIns, maxDel,maxSubRate, maxInDelRate)): # read was originally mapped but didn't pass the threshold
            self.flag += 0x4 # change its own "unmapped" flag
            if self.mate is not None:
                self.mate.flag += 0x8 # change mate's "mate_unmapped" flag
                self.cigarstring = '*'
        if self.is_paired() and self.mate is not None and not self.mate_is_unmapped() and self.mate.is_read_bad(maxSub, maxIns, maxDel,maxSubRate, maxInDelRate): # mate was originally mapped but didn't pass the threshold
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
        if self.is_paired() and self.mate is not None and  not self.mate_is_unmapped() and self.mate.is_read_bad(maxSub, maxIns, maxDel,maxSubRate, maxInDelRate): # mate was originally mapped but didn't pass the threshold
            self.fields[7] = self.mate.rstart
        if self.mate is not None:
            self.mate.fields[7] = self.rstart
        # [8] TLEN stays the same. signed observed Template LENgth.
        # If all segments are mapped to the same reference, the unsigned observed template length equals the number of bases from the leftmost mapped base to the rightmost mapped base. 
        # The leftmost segment has a plus sign and the rightmost has a minus sign. The sign of segments in the middle is undefined. It is set as 0 for single-segment template or when the information is unavailable. 
        # TODO: this could have changed
        # [9] SEQ segment SEQuence: trimmed original sequence (only the part mapped to PacBio sequence)
        self.trim_qseq()
        self.fields[5], self.fields[9] = self.cigarstring, self.qSeq
        if self.mate is not None:
            self.mate.trim_qseq()
            self.mate.fields[5], self.mate.fields[9] = self.mate.cigarstring, self.mate.qSeq
        # [10] QUAL stays the same

        ## Paste all fields together, separated by tabs, if at least one read of the pair is mapped. 
        ## Two lines if read has a mate, no matter if they are mapped or not
        if not self.is_unmapped() or (self.mate is not None and not self.mate_is_unmapped()):
            record = '\t'.join(map(str,self.fields))
            if self.mate is not None:
                record += '\n' + '\t'.join(self.mate.fields)
        else:
            record = ''
        return record

    # re-align read to PacBio sequence and shift indels to the leftmost possible positions
    def re_align(self, rseq, maxSub=-1, maxIns=-1, maxDel=-1,maxSubRate=0.05, maxInDelRate=0.3, max_round = 10):
        ''' realign the read to PacBio sequence, including retrieving clipped parts and scrubbing '''
        #print self.qname
        done = False # indicator of whether realigning is done
        rLen = len(rseq)
        while not done:
            #print "re-align"
            keep_left, keep_right,seg = self.get_seg(rseq) # get the correct segment for realigning

            ref_region_start = max( self.rstart - keep_left - 5, 1)
            ref_region_end = min(self.get_rend() + keep_right + 5, rLen)
            #print ref_region_start, ref_region_end
            #print "sequence A: ", rseq[(ref_region_start-1):ref_region_end]
            #print "sequence B: ", self.qSeq
            #print "\n"
            ## global pairwise alignment, 1 penalty for mismatch and 0.9 penalty for indels, opening and ending gaps in Illumina reads don't get penalized
            realign_res = pairwise2.align.globalms(rseq[(ref_region_start-1):ref_region_end], seg, 0, -1, -0.9, -0.9, penalize_end_gaps=[True, False])
            #print len(realign_res), " equivalent good mappings, after extending the Illumina read"
            #print format_alignment(*realign_res[0])
            new_align, align_start, align_end = metalrec_lib.pick_align(realign_res) # pick the best mapping: indel positions are the leftmost collectively
            #align_start = new_align[3] # starting position of the alignment for the new extended alignment
            #print format_alignment(*new_align) # for DEBUG
            new_align1 = metalrec_lib.shift_to_left_chop(new_align)
            rounds = 0
            while new_align1[:2] != new_align[:2] and rounds < max_round:
                #print "realign round ", rounds
                #print format_alignment(*new_align1) # for DEBUG
                new_align = list(new_align1)
                #print "before realign: ", new_align
                new_align1 = metalrec_lib.shift_to_left_chop(new_align)
                #print "after  realign: ", new_align1
                #print format_alignment(*new_align1)
                rounds += 1
                
            #print "done"
            #print format_alignment(*new_align) # for DEBUG

            # update information in the sam record
            self.rstart = ref_region_start + align_start # starting position
            new_cigarstring,first_non_gap,last_non_gap = metalrec_lib.get_cigar(new_align1[0], new_align1[1]) # get the cigar string for the new alignment
            cigar_info = metalrec_lib.cigar(self.cigarstring)
            left_clip_len = cigar_info['left_clip_len']
            right_clip_len = cigar_info['right_clip_len']
            # add the clipped part back to the cigar string
            if keep_left < left_clip_len:
                new_cigarstring = str(left_clip_len - keep_left) + 'S' + new_cigarstring
            if keep_right < right_clip_len:
                new_cigarstring = new_cigarstring + str(right_clip_len - keep_right) + 'S'
            self.cigarstring = new_cigarstring
            self.fields[3] = str(self.rstart)
            self.fields[5] = self.cigarstring
            self.alignRecord = '\t'.join(self.fields)
            #print new_cigarstring

            # check if more retrieving needs to be done
            if self.cigarstring.find('S') == -1:
                done = True
            else:
                keep_left, keep_right,seg = self.get_seg(rseq) # get the correct segment for realigning
                if keep_left == 0 and keep_right == 0:
                    done = True

            # check if there are already too many errors, if so, mark the read to be bad, and return empty dictionaries 
            if self.is_read_bad(maxSub, maxIns, maxDel,maxSubRate, maxInDelRate):
                #print "too many errors, discard"
                done = True
                return dict(), dict()
        #print "Done"        
        pos_dict, ins_dict = metalrec_lib.get_bases_from_align(new_align1, ref_region_start + align_start)

        return pos_dict, ins_dict
