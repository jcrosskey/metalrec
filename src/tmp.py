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
import metalrec_lib
rseq = metalrec_lib.read_single_seq(ref_fasta)
ref_bps, ref_ins_dict, read_info = metalrec_lib.read_and_process_sam(samfile, rseq, maxSubRate=0.1,outsam='a.sam')
good_regions = metalrec_lib.get_good_regions(ref_bps, rseq, minPacBioLen=1000, minCV=1)
poly_bps, poly_ins, consensus_bps, consensus_ins, cvs = metalrec_lib.get_poly_pos(ref_bps, ref_ins_dict, good_regions[0])
newSeq, bp_pos_dict, ins_pos_dict = metalrec_lib.ref_extension(poly_bps, poly_ins, consensus_bps, consensus_ins, rseq)
poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext = metalrec_lib.update_pos_info(poly_bps, poly_ins, consensus_bps, consensus_ins, bp_pos_dict, ins_pos_dict)
type_array = metalrec_lib.make_type_array(poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext)
ref_array = metalrec_lib.make_ref_array(consensus_bps_ext, consensus_ins_ext, type_array)
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
metalrec_lib.get_reads_for_gap(read_array, (gap_start_ind0[0], gap_end_ind0[0]), skip_reads=Cvec0) # reads to fill the gap

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


