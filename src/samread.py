#!/usr/bin/python

import metalrec_lib
import time
import re
import sys
sys.path.append("/lustre/atlas/scratch/chaij1/csc124/biopython-1.64/lib/python2.7/site-packages")
import math
from Bio import pairwise2 # pairwise alignment using dynamic programming
from Bio.pairwise2 import format_alignment
from numpy import *

''' reverse complement function '''
revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','-':'-'}.get(B,'N') for B in x ][::-1]) # find reverse complement of a DNA sequence


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
        #self.trim_qseq()
        self.fields[9] = self.qSeq
        if self.mate is not None:
            #self.mate.trim_qseq()
            self.mate.fields[9] = self.mate.qSeq
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
    def re_align(self, rseq, maxSub=-1, maxIns=-1, maxDel=-1,maxSubRate=0.05, maxInDelRate=0.3, max_round = 10,checkEnds=True):
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
            s_time = time.time()
            realign_res = pairwise2.align.globalms(rseq[(ref_region_start-1):ref_region_end], seg, 0, -1, -0.9, -0.9, penalize_end_gaps=[True, False])
            e_time = time.time()
            sys.stdout.write("align time :" + str(e_time - s_time) +  " seconds\n")
            
            #print len(realign_res), " equivalent good mappings, after extending the Illumina read"
            #print format_alignment(*realign_res[0])
            new_align, align_start, align_end = metalrec_lib.pick_align(realign_res) # pick the best mapping: indel positions are the leftmost collectively
            #align_start = new_align[3] # starting position of the alignment for the new extended alignment
            #print format_alignment(*new_align) # for DEBUG
            new_align1 = metalrec_lib.shift_to_left_chop(new_align)
            rounds = 0
            s_time = time.time()
            while new_align1[:2] != new_align[:2] and rounds < max_round:
                #print "realign round ", rounds
                #print format_alignment(*new_align1) # for DEBUG
                new_align = list(new_align1)
                #print "before realign: ", new_align
                new_align1 = metalrec_lib.shift_to_left_chop(new_align)
                #print "after  realign: ", new_align1
                #print format_alignment(*new_align1)
                rounds += 1
                
            e_time = time.time()
            sys.stdout.write("shuffle time :" + str(e_time - s_time) +  " seconds\n")
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
        if checkEnds:
            left_trim, right_trim, rstart_shift = self.check_ends(maxSubRate=maxSubRate)
            #sys.stdout.write("left_trim: {}; right_trim: {},  rstart_shift: {}\n\n".format(left_trim, right_trim, rstart_shift))
            seq1 = new_align1[0][left_trim:]
            seq2 = new_align1[1][left_trim:]
            if right_trim > 0:
                seq1 = seq1[:-right_trim]
                seq2 = seq2[:-right_trim]
            new_align1 = (seq1, seq2, 0, len(seq1), new_align1[-1])
            # change the starting position on ref sequence
            self.rstart = self.rstart + rstart_shift
            self.fields[3] = str(self.rstart)

        pos_dict, ins_dict = metalrec_lib.get_bases_from_align(new_align1, self.rstart)

        return pos_dict, ins_dict

    # check the ends of read and see if there are too many errors
    def check_ends(self, maxSubRate=0.05):
        cigar_info = metalrec_lib.cigar(self.cigarstring)
        left_trim_len = 0
        left_clip_len = 0
        right_trim_len = 0
        right_clip_len = 0
        rstart_shift = 0
        if cigar_info['sub_len'] > 1: # only do this if there is more than 1 substitution errors
            if maxSubRate is None:
                maxSubRate = cigar_info['sub_len'] / float(cigar_info['seq_len']) # average substitution rate in this read
            all_Nums = array(map(int,re.findall('\d+',self.cigarstring)))
            all_Chars = array(re.findall('\D+',self.cigarstring))
            non_clip_chars = where(all_Chars!='S')[0]
            Nums = all_Nums[non_clip_chars]
            Chars = all_Chars[non_clip_chars]
            sub_Inds = where(Chars=='X')[0] # index of substitutions in Chars
            if len(sub_Inds) > 1:
                # left end
                subs = cumsum(Nums[sub_Inds]) # cumulative number of substitutions
                bases = cumsum(Nums)[sub_Inds[1:] - 1] # number of bases corresponding to the number of substitutions

                subRate = true_divide(subs[:-1], bases) # substitution rates up to each substitution error place
                ManySubInds = where(subRate > maxSubRate*1.15)[0] # indices where subRate is too high
                if len(ManySubInds) > 0:
                    leftmost_ind = sub_Inds[amax(ManySubInds)]
                    left_trim_len = cumsum(Nums)[leftmost_ind]
                    left_clip_len = left_trim_len - sum(Nums[where(Chars[:leftmost_ind] =='D')[0]])
                    rstart_shift = left_trim_len - sum(Nums[where(Chars[:leftmost_ind] =='I')[0]])

                # right end
                rNums = Nums[::-1]
                rChars = Chars[::-1]
                rsub_Inds = where(rChars=='X')[0] # index of substitutions in Chars
                subs = cumsum(rNums[rsub_Inds]) # cumulative number of substitutions
                bases = cumsum(rNums)[rsub_Inds[1:] - 1] # number of bases corresponding to the number of substitutions

                subRate = true_divide(subs[:-1], bases) # substitution rates up to each substitution error place
                ManySubInds = where(subRate > maxSubRate*1.15)[0] # indices where subRate is too high
                if len(ManySubInds) > 0:
                    rightmost_ind = rsub_Inds[amax(ManySubInds)]
                    right_trim_len = cumsum(rNums)[ rightmost_ind ]
                    right_clip_len = right_trim_len - sum(rNums[where(rChars[:rightmost_ind ] =='D')[0]])
            #sys.stdout.write("left_clip_len: {} right_clip_len: {}\n".format(left_clip_len, right_clip_len))
            if left_trim_len > 0:
                Nums = Nums[(leftmost_ind + 1):]
                Chars = Chars[(leftmost_ind + 1):]
            if right_trim_len > 0:
                Nums = Nums[: -(rightmost_ind+1)]
                Chars = Chars[ : -(rightmost_ind+1)]
            Nums = Nums.tolist()
            Chars = Chars.tolist()
            if all_Chars[0] == 'S':
                Nums.insert(0, left_clip_len + all_Nums[0])
                Chars.insert(0, 'S')
            elif left_clip_len > 0:
                Nums.insert(0, left_clip_len)
                Chars.insert(0, 'S')

            if all_Chars[-1] == 'S':
                Nums.append(right_clip_len + all_Nums[-1])
                Chars.append('S')
            elif right_clip_len > 0:
                Nums.append(right_clip_len)
                Chars.append('S')

            res = [None] * (len(Nums) + len(Chars))
            res[::2] = map(str,Nums)
            res[1::2] = Chars
            self.cigarstring = ''.join(res)
            self.fields[5] = self.cigarstring
        return left_trim_len, right_trim_len, rstart_shift

class BlasrRead:
    ''' Class representing an aligned read, particularly from sam file, instead of bam file '''
    def __init__(self, BlasrRecord):
        self.BlasrRecord = BlasrRecord.strip('\n')
        fields = BlasrRecord.split()
        self.fields = fields
        self.qname = fields[0] # Query template name
        self.qLen = int(fields[1])
        self.qStart = int(fields[2])
        self.qEnd = int(fields[3])
        self.qStrand = fields[4]
        self.rName = fields[5]
        self.rLength = int(fields[6])
        self.rStart = int(fields[7])
        self.rEnd = int(fields[8])
        self.rStrand = fields[9]
        self.score = int(fields[10])
        self.nMatch = int(fields[11])
        self.nMisMatch = int(fields[12])
        self.nIns = int(fields[13])
        self.nDel = int(fields[14])
        self.mapQV = int(fields[15])
        self.qSeq = fields[16] if self.rStrand == '+' else revcompl(fields[16])
        self.rSeq = fields[18] if self.rStrand == '+' else revcompl(fields[18])
        self.mapLen = self.qEnd - self.qStart + self.nDel
        self.percentSimilarity = float(self.nMatch) / float(self.mapLen) # percent similarity


    # check and see if this mapping is too noisy to be included, if so, change the flag indicating if the read was mapped
    def is_read_bad(self, maxSub=-1, maxIns=-1, maxDel=-1,maxSubRate=0.05, maxInDelRate=0.3):
        indel_len = self.nIns + self.nDel
        # substitution rate is relative to the aligned region length on the PacBio sequence
        # indelRate is relative to the total length of the mapped segment of the Illumina sequence
        if self.nMisMatch > maxSubRate * self.mapLen or indel_len > maxInDelRate * self.mapLen: # subRate * mapped ref region; indelRate * seqLen 
            #print "sub_len: ", cigar_info['sub_len'], ", indel len: ", indel_len #DEBUG
            return True
        # Finally, if it passes all the thresholds, it's a good record
        else:
            return False

    # re-align read to PacBio sequence and shift indels to the leftmost possible positions
    def re_align(self, maxSub=-1, maxIns=-1, maxDel=-1,maxSubRate=0.05, maxInDelRate=0.3, max_round = 10,checkEnds=True):
        ''' realign the read to PacBio sequence, including retrieving clipped parts and scrubbing '''
        new_align = (self.rSeq, self.qSeq, 0, 0, len(self.qSeq)) # same format as output of pairwise... tuple of 5 entries
        #print format_alignment(*new_align) # for DEBUG
        new_align1 = metalrec_lib.shift_to_left_chop(new_align)
        rounds = 0
        #s_time = time.time()
        while new_align1[:2] != new_align[:2] and rounds < max_round:
            #print "realign round ", rounds
            #print format_alignment(*new_align1) # for DEBUG
            new_align = list(new_align1)
            #print "before realign: ", new_align
            new_align1 = metalrec_lib.shift_to_left_chop(new_align)
            #print "after  realign: ", new_align1
            #print format_alignment(*new_align1)
            rounds += 1
                
        #e_time = time.time()
        #sys.stdout.write("shuffle time :" + str(e_time - s_time) +  " seconds\n")
        #print "done"
        #print format_alignment(*new_align) # for DEBUG
        self.qSeq = new_align1[1]
        self.rSeq = new_align1[0]

        # check if there are already too many errors, if so, mark the read to be bad, and return empty dictionaries 
        if self.is_read_bad(maxSub, maxIns, maxDel,maxSubRate, maxInDelRate):
            #print "too many errors, discard"
            return dict(), dict()
        else:
            pos_dict, ins_dict = metalrec_lib.get_bases_from_align(new_align1, self.rStart + 1)
            return pos_dict, ins_dict


    def generate_sam_record(self, maxSub=-1, maxIns=-1, maxDel=-1,maxSubRate=0.1, maxInDelRate=0.3):
        # [0] qname stays the same
        # [1] new flag: mapped/unmapped status might change
        self.flag = 0
        if self.qStrand != self.rStrand:
            self.flag += 0x8 # reverse complement, other flag bits are not set for now (only for single end reads)
        # [2] rname stays the same
        # [3] pos: 1-based leftmost mapping POSition might change after re-mapping 
        # [4] MAPQ stays the same
        # [5] CIGAR string might have changed after re-mapping to PacBio sequence
        self.cigar, first_non_gap, last_non_gap = metalrec_lib.get_cigar(self.rSeq, self.qSeq)
        #print self.cigar, first_non_gap, last_non_gap
        if first_non_gap != 0:
            self.rStart += first_non_gap
        #if self.qStart != 0:
        #    self.cigar = str(self.qStart) + 'S' + self.cigar
        #if self.qEnd != self.qLen:
        #    self.cigar = self.cigar + str(self.qLen - self.qEnd) + 'S'
        # [6] RNEXT Ref. name of the mate/next read: should always be = since there is only 1 PacBio sequence 
        # [7] PNEXT Position of the mate/next read: might change because of re-mapping
            # will be set to 0 now
        # [8] TLEN stays the same. signed observed Template LENgth. set to 0 for now
        # If all segments are mapped to the same reference, the unsigned observed template length equals the number of bases from the leftmost mapped base to the rightmost mapped base. 
        # The leftmost segment has a plus sign and the rightmost has a minus sign. The sign of segments in the middle is undefined. It is set as 0 for single-segment template or when the information is unavailable. 
        # [9] SEQ segment SEQuence: trimmed original sequence (only the part mapped to PacBio sequence)
        # [10] QUAL stays the same

        ## Paste all fields together, separated by tabs, if at least one read of the pair is mapped. 
        record = '\t'.join([self.qname, str(self.flag), self.rName, str(self.rStart + 1), str(self.mapQV), self.cigar, "=", "0", "0", re.sub('-', '',self.qSeq),'*']) + '\n'
        return record
