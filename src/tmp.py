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
