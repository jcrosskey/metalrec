# -*- coding: utf-8 -*-
"""
Created on Fri Jun 27 14:08:08 2014

@author: cjg
"""

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
rseq = 'CCGTTCGGCGTCAAAGCCCACGCGCAGCTCTTCAGCGGTCCGCGTGGCGCGGGCGCAGCATGGGTGGCGGTCGCTACGGAAGGCACTACTAGCTACTAGCGGCGAGAGCGAGGGTAGGCGGGGCGAGCTGGCAGGTTGGCGGGTACGCGCGCCGCCATTGCAATACGCTTGCCGTCGCGGCGGGCCTTCGAGCGTGGCATG'
qseq = 'TGGGCGCGGGCGCAGCATGGGTGGCGGTCGCTACGGAAGGCACTACTAGCTACTAGCGGCGAGAGCGAGGGTAGGCGGGCGAGCTGGCAGGTTGGCGGGTACGCCGCGCCCGCCATTGCAATACGCTTGCCG'
a = pairwise2.align.globalms(rseq, qseq, 0, -1, -0.9,-0.9, penalize_end_gaps=[True, False])

for i in a:
    print format_alignment(*i)
    
#samfile = "/Users/cjg/Work/PacBio/Results/MockCommunity/Illu_reads_to_roi_each/m130828_015813_00123_c100564312550000001823090912221380_s1_p0__151312__ccs/bbmap/test.sam"
#ref_fasta = "/Users/cjg/Work/PacBio/Data/MockCommunity/PacBio/RoI_90_single/m130828_015813_00123_c100564312550000001823090912221380_s1_p0__151312__ccs.fasta"
samfile = "/Users/cjg/Work/PacBio/Results/Wetlands/subreads_with_good_regions/m131016_052225_00123_c100575992550000001823095504021421_s1_p0__54538__13380_14430/a.sam"
ref_fasta = "/Users/cjg/Work/PacBio/Results/Wetlands/subreads_with_good_regions/m131016_052225_00123_c100575992550000001823095504021421_s1_p0__54538__13380_14430/m131016_052225_00123_c100575992550000001823095504021421_s1_p0__54538__13380_14430.fasta"

samfile = "/Users/cjg/Work/PacBio/metalrec/test/bbmap_red.sam"
ref_fasta = "/Users/cjg/Work/PacBio/metalrec/test/m130828_041445_00123_c100564312550000001823090912221381_s1_p0__54536__2157_4456.fasta"
import metalrec_lib
rseq = metalrec_lib.read_single_seq(ref_fasta)
# process sam file and save the read info, considering each read as single end read
ref_bps, ref_ins_dict, read_info = metalrec_lib.read_and_process_sam(samfile, rseq, maxSubRate=0.1,outsam='/Users/cjg/Desktop/m131016_052225_00123_c100575992550000001823095504021421_s1_p0__54542__11111_12832/a.sam')
# process sam file and save the read info, considering pair-end reads
ref_bps, ref_ins_dict, read_info = metalrec_lib.read_and_process_sam_pair(samfile, rseq, maxSubRate=0.1,outsam='/Users/cjg/Desktop/m131016_052225_00123_c100575992550000001823095504021421_s1_p0__54542__11111_12832/a_p.sam')

# key for the pair of mates with ID HISEQ11:283:H97Y1ADXX:2:2205:12196:34085 in read_info_p
for key in read_info_p:
    if key[:15] == '42G43C44A45A46C':
        print key, '\n'
        a = key


a = a.split(':') # separate the non-insertion and insertion positions

good_regions = metalrec_lib.get_good_regions(ref_bps, rseq, minPacBioLen=1000, minCV=1)
poly_bps, poly_ins, consensus_bps, consensus_ins, cvs = metalrec_lib.get_poly_pos(ref_bps, ref_ins_dict, good_regions[0])
newSeq, bp_pos_dict, ins_pos_dict = metalrec_lib.ref_extension(poly_bps, poly_ins, consensus_bps, consensus_ins, rseq)
poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext = metalrec_lib.update_pos_info(poly_bps, poly_ins, consensus_bps, consensus_ins, bp_pos_dict, ins_pos_dict)
type_array = metalrec_lib.make_type_array(poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext)
#ref_array = metalrec_lib.make_ref_array(consensus_bps_ext, consensus_ins_ext, type_array)
read_array, read_counts = metalrec_lib.make_read_array(read_info, bp_pos_dict, ins_pos_dict, type_array, poly_bps, poly_ins, consensus_bps, consensus_ins)

a = read_info.keys()[1]
array_a = metalrec_lib.make_read_array1d(a, bp_pos_dict, ins_pos_dict, type_array,poly_bps, poly_ins, consensus_bps, consensus_ins)
array_b = metalrec_lib.make_read_array1d(read_info.keys()[0], bp_pos_dict, ins_pos_dict, type_array,poly_bps, poly_ins, consensus_bps, consensus_ins)
aa = array_a.reshape(-1,5)
where(array_a[4::5]==1)[0] # positions where there is a 'D'
hstack((where(type_array==3)[0], where(type_array==1)[0])) # insertion positions

adjac = gaps[1:] - gaps[:-1] # difference between a position in the gap vec and the previous position, if the difference is 1, then it's consecutive gap
gap_starts = where(adjac != 1)[0] # positions where new gap starts
gap_lens = gap_starts -  concatenate(( array([0]), gap_starts[:-1] )) # length of the gaps ( not including the last gap)
gap_lens = concatenate( (gap_lens, array([len(gaps) - gap_starts[-1] - 1]) )) # append the length of the last gap
max_gap_ind = argmax(gap_lens) # index of the widest gap

ref0 = metalrec_lib.get_consensus_from_array(read_array)
Cvec0 = metalrec_lib.get_compatible_reads(ref0, read_array)
gap_pos0 = metalrec_lib.gap_pos(ref0, read_array, Cvec0)
gap_start_ind0, gap_end_ind0 = metalrec_lib.get_gaps(gap_pos0)
reads_ind, reads_cov = metalrec_lib.get_reads_for_gap(read_array, (gap_start_ind0[0], gap_end_ind0[0]), skip_reads=Cvec0) # reads to fill the gap, and how many bases they can fill

ref1 = metalrec_lib.get_new_ref(ref0, 10, read_array) # get a new ref
Cvec1 = metalrec_lib.get_compatible_reads(ref1, read_array)
gap_pos1 = metalrec_lib.gap_pos(ref1, read_array, Cvec1)
gap_start_ind1, gap_end_ind1 = metalrec_lib.get_gaps(gap_pos1)
metalrec_lib.get_reads_for_gap(read_array, (gap_start_ind1[0], gap_end_ind1[0]), skip_reads=Cvec1) # reads to fill the gap

ref2 = metalrec_lib.get_new_ref(ref0, 319, read_array) # get a new ref
Cvec2 = metalrec_lib.get_compatible_reads(ref2, read_array)
gap_pos2 = metalrec_lib.gap_pos(ref2, read_array, Cvec2)
gap_start_ind2, gap_end_ind2 = metalrec_lib.get_gaps(gap_pos2)
metalrec_lib.get_reads_for_gap(read_array, (gap_start_ind2[0], gap_end_ind2[0]), skip_reads=Cvec2) # reads to fill the gap


# check for the position in the existing PacBio sequence
for k, z in bp_pos_dict.items():
    if z == 251:
        print k
        break

## test of newly added functions in samread class ##
import samread
import metalrec_lib
from numpy import *
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

samfile = "/Users/cjg/Work/PacBio/Results/MockCommunity/02_Debug/m130828_015813_00123_c100564312550000001823090912221380_s1_p0__117926__532_2975/m130828_015813_00123_c100564312550000001823090912221380_s1_p0__117926__532_2975.sam"
ref_fasta = "/Users/cjg/Work/PacBio/Results/MockCommunity/02_Debug/m130828_015813_00123_c100564312550000001823090912221380_s1_p0__117926__532_2975/m130828_015813_00123_c100564312550000001823090912221380_s1_p0__117926__532_2975.fasta"
rseq = metalrec_lib.read_single_seq(ref_fasta)

reload(metalrec_lib)
reload(samread)
r = samread.SamRead(a)
r.re_align(rseq)

samIn = open(samfile,'r')
for i in xrange(4):
    a = samIn.readline()

while a.split('\t')[0] != 'HISEQ11:285:H987LADXX:2:2104:15628:91916':
    a = samIn.readline()


r = samread.SamRead(a) # SamRead object
samIn.close()

r1 = samread.SamRead(a) # SamRead object
rLen = len(rseq)
ref_region_start = max( r1.rstart - 5, 1)
ref_region_end = min(r1.get_rend() + 5, rLen)
#print r1.qSeq
r1.trim_qseq()
#print r1.qSeq
realign_res = pairwise2.align.globalms(rseq[(ref_region_start-1):ref_region_end], r1.qSeq, 0, -1, -0.9, -0.9, penalize_end_gaps=[True, False])
new_align = metalrec_lib.pick_align(realign_res) # pick the best mapping: indel positions are the leftmost collectively
#print format_alignment(*new_align) # for DEBUG
pos_dict, ins_dict = metalrec_lib.get_bases_from_align(new_align, ref_region_start + new_align[3])

ref_bps, ref_ins_dict, readinfo = metalrec_lib.read_and_process_sam_samread(samfile, rseq, maxSub=-1, maxDel=-1, maxIns=-1, maxSubRate=0.05, maxInDelRate=0.3, verbose=True)
good_regions, cov_bps, avg_cov_depth= metalrec_lib.get_good_regions(ref_bps, rseq, minGoodLen=400, minCV=1)
poly_bps, poly_ins, consensus_bps, consensus_ins, cvs = metalrec_lib.get_poly_pos(ref_bps, ref_ins_dict, good_regions[0])
newSeq, bp_pos_dict, ins_pos_dict = metalrec_lib.ref_extension(poly_bps, poly_ins, consensus_bps, consensus_ins, rseq, print_width=70, verbose=True)
poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext = metalrec_lib.update_pos_info(poly_bps, poly_ins, consensus_bps, consensus_ins, bp_pos_dict, ins_pos_dict)
type_array, ext_region = metalrec_lib.make_type_array(poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext,verbose=True)
#consensus_array = metalrec_lib.make_ref_array(consensus_bps_ext, consensus_ins_ext, type_array,ext_region)
ar = metalrec_lib.make_read_array1d(readinfo.keys()[0], bp_pos_dict, ins_pos_dict, type_array, poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext) 
short_r, long_r = metalrec_lib.array_to_seq(ar) # get the read sequence from its array
read_array, read_counts = metalrec_lib.make_read_array(readinfo, bp_pos_dict, ins_pos_dict, type_array, poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext)
ref0, max_len,gap_pos, Cvec = metalrec_lib.greedy_fill_gap(read_array, ref0=None, minOverlap=10,minOverlapRatio=0.1, verbose=True)
overlap_mat = metalrec_lib.get_overlapLen(ref0, read_array, Cvec=Cvec)
ref_new = metalrec_lib.fill_gap(read_array, 15, "/Users/cjg/Work/PacBio/Results/MockCommunity/02_Debug/m130828_015813_00123_c100564312550000001823090912221380_s1_p0__100124__481_2720/EC/goodreads.fasta", "/Users/cjg/Work/PacBio/Results/MockCommunity/02_Debug/m130828_015813_00123_c100564312550000001823090912221380_s1_p0__100124__481_2720/EC/region0", readinfo, verbose=True)
ref_new = metalrec_lib.fill_gap(read_array, "/Users/cjg/Work/PacBio/metalrec/test/good_reads.fasta", "/Users/cjg/Work/PacBio/metalrec/test", readinfo, verbose=False)
ref_new_short, ref_new_long = metalrec_lib.array_to_seq(ref_new)
new_ref = metalrec_lib.fill_gap(read_array, outFastaFile="/Users/cjg/Work/PacBio/metalrec/test/good_reads.fasta", outDir=None, readinfo=readinfo, verbose=True)
