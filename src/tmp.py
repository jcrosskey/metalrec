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
    
samfile = "/Users/cjg/Work/PacBio/Results/MockCommunity/Illu_reads_to_roi_each/m130828_015813_00123_c100564312550000001823090912221380_s1_p0__151312__ccs/bbmap/test.sam"
ref_fasta = "/Users/cjg/Work/PacBio/Data/MockCommunity/PacBio/RoI_90_single/m130828_015813_00123_c100564312550000001823090912221380_s1_p0__151312__ccs.fasta"
import metalrec_lib
rseq = metalrec_lib.read_single_seq(ref_fasta)
ref_bps, ref_ins_dict = metalrec_lib.read_and_process_sam(samfile, rseq)
good_regions = metalrec_lib.get_good_regions(ref_bps, ref_ins_dict, rseq, minPacBioLen=1000, minCV=10)
poly_bps, poly_ins, consensus_bps, consensus_ins = metalrec_lib.get_poly_pos(ref_bps, ref_ins_dict, good_regions[0])
ref_bps, ref_ins_dict, readinfo = metalrec_lib.read_and_process_sam(samfile, rseq)
