#!/usr/bin/python

''' Collection of functions used in MetaLREC project
'''
import re
import sys
import samread
from Bio import pairwise2 # pairwise alignment using dynamic programming
from Bio.pairwise2 import format_alignment

## ======================================================================
## From CIGAR string, find number of indels and substitutions  
## also the length of the read, the sequence, and the mapped reference region 
## arguments: CIGAR string
## ======================================================================
def cigar(cigar_string):
    """ parse CIGAR string from .sam file, find number of matches, mismatches, insertion (to ref), deletion,
        soft clipping, hard clipping, padding
        return: 1-level dictionary
    """
    char = re.findall('\D',cigar_string) # operation characters, MIDNSHP=X
    char = [x.upper() for x in char]  # convert to upper case
    count = map(int,re.findall('\d+',cigar_string))
    read_len = sum(count[i] for i in xrange(len(char)) if char[i] in 'MIS=XH') # length of the mapped read (including clipped part)
    seq_len = sum(count[i] for i in xrange(len(char)) if char[i] in 'MI=X') # length of the mapped read (excluding clipped part)
    ref_region_len = sum(count[i] for i in xrange(len(char)) if char[i] in 'MD=X') # length of the region on the reference sequence corresponding to the read
    ins_len = sum(count[i] for i in xrange(len(char)) if char[i] == 'I') # insert bps
    del_len = sum(count[i] for i in xrange(len(char)) if char[i] == 'D') # deletion bps
    try:
        max_ins_len = max(count[i] for i in xrange(len(char)) if char[i] == 'I') #  maximum insertion stretch
    except ValueError:
        max_ins_len = 0

    try:
        max_del_len = max(count[i] for i in xrange(len(char)) if char[i] == 'D') #  maximum deletion stretch
    except ValueError:
        max_del_len = 0

    match_len = None # these are only available with X= in the cigar string, sam 1.4 and after
    sub_len = None
    max_sub_len = 0
    if '=' in char or 'X' in char: # if either '=' or 'X' is detected, indicating sam 1.4 format for CIGAR string
        match_len = sum(count[i] for i in xrange(len(char)) if char[i] == '=') # matched bps
        sub_len = sum(count[i] for i in xrange(len(char)) if char[i] == 'X') # substitution bps
        try:
            max_sub_len = max(count[i] for i in xrange(len(char)) if char[i] == 'X') # maximum substitution stretch
        except ValueError:
            max_sub_len = 0
        align_len = match_len + sub_len
    else:   # sam 1.3, aligned bps from 'M'
        align_len = sum(count[i] for i in xrange(len(char)) if char[i] == 'M') # aligned bps(including matching and substitution)

    pad_len = sum(count[i] for i in xrange(len(char)) if char[i] == 'P') # padded bps, inserted both in reference sequence and the read sequence)
    return {'read_len':read_len, 'seq_len': seq_len, 'ref_len':ref_region_len, 'ins_len':ins_len, 'del_len':del_len, 'match_len':match_len, 'sub_len':sub_len, 'align_len':align_len, 'pad_len':pad_len, 'max_ins': max_ins_len, 'max_del': max_del_len, 'max_sub': max_sub_len}

## ======================================================================
## From MD tag, find the number of matched, deletion, and substitution bps 
## arguments: MD tag
## ======================================================================
def md(MD_tag):
    """
    Given MD tag and a sequence, find the number of matched, deletion, and substitution bps 

    Return: 1-level dictionary 
    """
    #matchLens = re.sub("\D"," ",MD_tag).split() #replace non-digit chars by blanks
    #matchLens = map(int, matchLens) #length of the matches before a mismatch/deletion
    match_len = sum(map(int,re.findall("\d+",MD_tag))) #Get all the numbers in MD tag, these are the number of matches
    nonMatches = re.findall("\D+",MD_tag) #mismatch/deletion, replace digits by blanks
    del_len = 0
    sub_len = 0
    for string in nonMatches:
        if string.startswith('^'): # deletion
            del_len += (len(string) - 1)
        else: # substitution
            sub_len += len(string)
    return {'match_len':match_len, 'del_len':del_len, 'sub_len':sub_len}

## ======================================================================
## function to check the CIGAR string to see if this record is bad or not
## works with sam format 1.4
## later should make it work with sam format 1.3
## ======================================================================
def is_record_bad(alignRecord,maxSub=3, maxIns=3, maxDel=3,maxSubRate=0.02, maxInsRate=0.2, maxDelRate=0.2):
    ''' Test and see if an alignment record is bad in the sam file.
        Input:  alignRecord - a mapping line in sam file
                maxSub - maximum stretches of substitution
                maxIns - maximum stretches of insertion
                maxDel - maximum stretches of deletion
                maxSubRate - maximum total substitution rate of the read
                maxInsRate - maximum total insertion rate of the read
                maxDelRate - maximum total deletion rate of the read
        Output: boolean value, True if bad else False
    '''
    fields = alignRecord.split("\t") # split by tabs
    cigarstring = fields[5] # CIGAR string
    if cigarstring == '*': # if cigar string is not available, treat as bad (no longer check the bitwise flag to make sure it's really not mapped)
        return True
    else:
        cigar_info = cigar(cigarstring)
        # if consecutive sub or indels is longer than the threshold, treat as bad
        if cigar_info['max_ins'] > maxIns or cigar_info['max_sub'] > maxSub or cigar_info['max_del'] > maxDel:
            return True
        # if any kind of error count exceeds the query sequence length * maximum allowed rate, also bad
        else:
            if cigar_info['sub_len'] > maxSubRate * cigar_info['seq_len']:
                return True
            elif cigar_info['ins_len'] > maxInsRate * cigar_info['seq_len']:
                return True
            elif cigar_info['del_len'] > maxDelRate * cigar_info['seq_len']:
                return True
            # Finally, if it passes all the thresholds, it's a good record
            else:
                return False
    
## ======================================================================
## Remove low quality alignments from short reads to long read, from the 
## mapping results.
## ======================================================================
def clean_samfile(samFile,samNew, rseq, maxSub=3, maxIns=3, maxDel=3,maxSubRate=0.02, maxInsRate=0.2, maxDelRate=0.2):
    newsam = open(samNew,'w')
    keepRec = 0
    discardRec = 0
    with open(samFile,'r') as mysam:
        for line in mysam:
            if line[0] == '@': # copy header lines
                newsam.write(line)
                if line[1:3] == 'SQ':
                    rLen = int(line[(line.find('LN:') + len('LN:')) : line.find('\t',line.find('LN:'))]) # reference sequence length
            else:
                record = line.strip('\n')
                fields = record.split('\t')
                if not is_record_bad(record, maxSub, maxIns, maxDel, maxSubRate, maxInsRate, maxDelRate): # if this alignment is good
                    myread = samread.SamRead(record)
                    # improve the alignment to the reference sequence by dynamic programming
                    # original starting and ending positions on the reference sequence of the mapped region
                    ref_region_start = max( myread.rstart - 10, 1)
                    ref_region_end = min(myread.get_rend() + 10, rLen)
                    # query sequence with clipped part trimmed
                    trimmed_qseq = myread.get_trim_qseq()
                    # redo global alignment using dynamic programming
                    realign_res = pairwise2.align.globalms(rseq[(ref_region_start-1):ref_region_end], trimmed_qseq, 0, -1, -0.9, -0.9, penalize_end_gaps=[True, False])
                    new_align = pick_align(realign_res) # pick the first mapping 
                    cigarstring,first_non_gap = get_cigar(new_align[0], new_align[1]) # get the cigar string for the new alignment
                    # update information in the sam record
                    fields[3] = str(ref_region_start + new_align[3]) # starting position
                    fields[5] = cigarstring # cigar string
                    fields[9] = trimmed_qseq # trimmed query sequence
                    # write updated record in the new file
                    line = '\t'.join(fields) + '\n'
                    newsam.write(line)
                    keepRec += 1
                else:
                    discardRec += 1
    newsam.close()
    # report some summary
    sys.stdout.write('Total number of records kept is {}. \n'.format(keepRec))
    sys.stdout.write('Total number of records discarded is {}. \n'.format(discardRec))

## ======================================================================
## From the CIGAR string of an alignment, the start position on the reference,
## and the mapped segment, find nucleotides base by base
## ======================================================================
def get_bases(cigar_string, qseq='', start_pos=''):
    ''' from CIGAR string, query segment (in alignment record), and starting position (1-based) on the ref sequence, return position wise base call from the read.
        Assuming that cigar_string is available, doesn't matter if it's sam 1.3 format or 1.4 format
        If there is only one input argument, treat it as the whole alignment record.
        
        Input: cigar_string, query aligned segment, and the 1-based starting mapping position on the reference sequence
        
        Output: (pos_dict, ins_dict), tuple of 2 dictionaries, one for the non-insertion positions, and one for the insertion positions.
                When there is a deletion from the reference sequence, the base called will be "D"
                pos_dict: ref_pos => (query_pos, base), or ref_pos => 'D'
                ins_dict: ref_pos => (query_pos, baese(s)) # the inserted length could be 1, or greater than 1
    '''
    # the whole alignment record
    if qseq == '' and start_pos == '':
        fields = cigar_string.strip('\n').split('\t')
        cigar_string = fields[5]
        qseq = fields[9]
        start_pos = int(fields[3])
    char = re.findall('\D',cigar_string) # operation characters, MIDNSHP=X
    char = [x.upper() for x in char]  # convert to upper case
    count = map(int,re.findall('\d+',cigar_string)) # corresponding count for each operation

    pos_dict = dict() # ref_pos (1-based) => base call
    ins_dict = dict() # ref_pos (1-based) => inserted base
    query_pos = 0

    # left end of the alignment
    if char[0] == 'S' or char[0] == 'N': # skipped region on the reference sequence (clipped or skipped)
        query_pos += count[0] # mapping starting position on the query segment (0-based)
        char = char[1:] # delete the first CIGAR operation and the corresponding count
        count = count[1:]
    elif char[0] == 'H': # hard clipping at the beginning
        char = char[1:]
        count = count[1:]

    # clipping (soft or hard) or skipping at the end of the mapping
    if char[-1] in 'SNH': 
        char = char[:-1]
        count = count[:-1]

    # Base calling, loop all the positions except clipping and skipping
    for pos in xrange(len(char)):
        if char[pos] == 'M' or char[pos] == 'X' or char[pos] == '=': # matching or mismatching
            for i in xrange(count[pos]):
                pos_dict[start_pos] = (query_pos, qseq[query_pos]) # base at the ref position in the read
                start_pos += 1
                query_pos += 1
        elif char[pos] == 'I': # insertion into reference sequence
            ins_dict[start_pos] = (query_pos, qseq[query_pos:(query_pos + count[pos])])
            query_pos += count[pos]
        elif char[pos] == 'D' or char[pos] == 'N': # deletion or skipped region from reference sequence
            for i in xrange(count[pos]):
                pos_dict[start_pos] = 'D'
                start_pos += 1
        else: # unknown CIGAR operations, print message
            sys.stdout.write(" unknown CIGAR operation: {}\n".format(char[pos]))

    return pos_dict, ins_dict

def get_bases_from_align(align, start_pos):
    ''' same as above, but from pairwise alignment results (of Bio.pairwise2 module), summarize the mapping '''
    pos_dict = dict() # ref_pos (1-based) => base call
    ins_dict = dict() # ref_pos (1-based) => inserted base
    
    seqA, seqB, score, begin, end = align
    ref_pos = start_pos
    query_pos = 0
    for i in xrange(len(seqA)):
        if seqA[i] == '-':
            ins_dict[ ref_pos ] = (query_pos, seqB[query_pos])
            query_pos += 1
        else:
            if seqB[i] == '-':
                pos_dict[ ref_pos ] = (query_pos, 'D')
            else:
                pos_dict[ ref_pos ] = (query_pos, seqB[query_pos])
                query_pos += 1
            ref_pos += 1

    return pos_dict, ins_dict
## ======================================================================
## From the mapping sam file and the reference sequence,
## get the mapping information for the reference sequence, base by base
## ======================================================================
def read_sam(samFile, maxSub=3, maxIns=3, maxDel=3,maxSubRate=0.02, maxInsRate=0.2, maxDelRate=0.2):
    ''' Get consensus sequence from alignments of short reads to a long read.

        Input: samFile - .sam file generated by mapping

        Output: ref_bps - list of lists, one for each position on the reference sequence (without padding)
                ref_ins_dict - dictionary of insertions, one for each insertion position found in the alignments.
    '''
    alphabet = 'ACGTD' # all the possible base pairs at a position
    ref_ins_dict = dict() # global insertion dictionary for the reference sequence
    keepRec = 0
    
    with open(samFile, 'r') as mysam:
        for line in mysam:
            if line[0] == '@': # header line
                if line[1:3] == 'SQ': # reference sequence dictionary
                    rname = line[(line.find('SN:') + len('SN:')) : line.find('\t',line.find('SN:'))] # reference sequence name
                    rLen = int(line[(line.find('LN:') + len('LN:')) : line.find('\t',line.find('LN:'))]) # reference sequence length
                    print "length of the reference is", rLen
                    ref_bps = [ [0] * 5 for x in xrange(rLen) ]  # list of lists, one list corresponding to each position on the reference sequence
            else:
                line = line.strip()
                if not is_record_bad(line, maxSub, maxIns, maxDel, maxSubRate, maxInsRate, maxDelRate): # if this alignment is good
                    keepRec += 1
                    fields = line.split("\t") # split by tabs
                    cigarstring = fields[5] # CIGAR string
                    qname = fields[0]
                    qseq = fields[9]
                    start_pos = int(fields[3]) # starting mapping position on the reference sequence, 1-based
                    pos_dict, ins_dict = get_bases(cigarstring, qseq, start_pos) # base calling from this read

                    for pos in pos_dict: # all the matching/mismatching/deletion positions
                        #print pos, '=>', pos_dict[pos]
                        ref_bps[pos-1][alphabet.find(pos_dict[pos][-1])] += 1 # update the corresponding base pair frequencies in the reference sequence

                    for ins in ins_dict: # all the insertion positions
                        if not ref_ins_dict.has_key(ins): # if this position has not appeared in the big insertion dictionary, initialize it
                            ref_ins_dict[ins] = [0] * 4 # length is 4 because there is no 'D' at insertion position
                        ref_ins_dict[ins][alphabet.find(ins_dict[ins][-1])] += 1 # if this position was already seen, just update frequencies of bases

                    if keepRec % 10000 == 0:
                        sys.stdout.write('  processed {} good records\n'.format(keepRec))
    
    return ref_bps, ref_ins_dict

def read_and_process_sam(samFile,rseq, maxSub=3, maxIns=3, maxDel=3,maxSubRate=0.02, maxInsRate=0.2, maxDelRate=0.2, minPacBioLen=1000, minCV=10):
    ''' Get consensus sequence from alignments of short reads to a long read, in the process, filter out bad reads and improve mapping

        Input:  samFile - .sam file generated by mapping
                threshold parameters

        Output: ref_bps - list of lists, one for each position on the reference sequence (without padding)
                ref_ins_dict - dictionary of insertions, one for each insertion position found in the alignments.
    '''
    alphabet = 'ACGTD' # all the possible base pairs at a position
    ref_ins_dict = dict() # global insertion dictionary for the reference sequence
    keepRec = 0
    lineNum = 0
    
    with open(samFile, 'r') as mysam:
        for line in mysam:
            lineNum += 1 
            if line[0] == '@': # header line
                if line[1:3] == 'SQ': # reference sequence dictionary
                    rname = line[(line.find('SN:') + len('SN:')) : line.find('\t',line.find('SN:'))] # reference sequence name
                    rLen = int(line[(line.find('LN:') + len('LN:')) : line.find('\t',line.find('LN:'))]) # reference sequence length
                    print "length of the reference is", rLen
                    if rLen < minPacBioLen:
                        return 0 
                    ref_bps = [ [0] * 5 for x in xrange(rLen) ]  # list of lists, one list corresponding to each position on the reference sequence
            else:
                record = line.strip('\n')
                fields = record.split('\t')
                if not is_record_bad(record, maxSub, maxIns, maxDel, maxSubRate, maxInsRate, maxDelRate): # if this alignment is good
                    keepRec += 1
                    myread = samread.SamRead(record)
                    # improve the alignment to the reference sequence by dynamic programming
                    # extend original starting and ending positions on the reference sequence of the mapped region by 5 bps on each end
                    ref_region_start = max( myread.rstart - 5, 1)
                    ref_region_end = min(myread.get_rend() + 5, rLen)
                    # query sequence with clipped part trimmed
                    trimmed_qseq = myread.get_trim_qseq()
                    # redo global alignment using dynamic programming
                    realign_res = pairwise2.align.globalms(rseq[(ref_region_start-1):ref_region_end], trimmed_qseq, 0, -1, -0.9, -0.9, penalize_end_gaps=[True, False])
                    new_align = pick_align(realign_res) # pick the first mapping 
                    pos_dict, ins_dict = get_bases_from_align(new_align, ref_region_start + new_align[3])
                    if myread.qname == 'HISEQ03:379:C2WP8ACXX:7:1203:19706:9668/1':
                        print new_align
                        print pos_dict
                        print ins_dict
                        print get_cigar(new_align[0], new_align[1])
                    for pos in pos_dict: # all the matching/mismatching/deletion positions
                        ref_bps[pos-1][alphabet.find(pos_dict[pos][-1])] += 1 # update the corresponding base pair frequencies in the reference sequence

                    for ins in ins_dict: # all the insertion positions
                        if ins not in ref_ins_dict: # if this position has not appeared in the big insertion dictionary, initialize it
                            ref_ins_dict[ins] = [0] * 4 # length is 4 because there is no 'D' at insertion position
                        ref_ins_dict[ins][alphabet.find(ins_dict[ins][-1])] += 1 # if this position was already seen, just update frequencies of bases

                    if keepRec % 1000 == 0:
                        sys.stdout.write('  processed {} good records\n'.format(keepRec))
                    if keepRec % 5000 == 0:
                        print lineNum
                        return ref_bps, ref_ins_dict
    
    return ref_bps, ref_ins_dict

## ======================================================================
## From the output of function get_consensus ( list of mat/mismat/del positions, and dictionary of ins positions),
## generate the consensus reference sequence and statistics.
## ======================================================================
def ref_consensus(ref_bps, ref_ins_dict, rSeq):
    ''' From the output of function get_consensus ( list of mat/mismat/del positions, and dictionary of ins positions), generate the consensus reference sequence and statistics.
        
        Input: ref_bps - list of lists, one for each position on the reference sequence (without padding)
               ref_ins_dict - dictionary of insertions, one for each insertion position found in the alignments.

        Output: A consensus sequence, and statistics...
    '''
    alphabet = 'ACGTD' # all the possible base pairs at a position
    rLen = len(rSeq)
    cov_depths = []
    true_ins = dict() # 'true' insertions into reference sequence
    # loop through all the positions of the reference sequence
    for i in xrange(len(rSeq)):
        non_ins_bps = ref_bps[i] # non-insertion base calls
        cov_depth = sum(non_ins_bps) # coverage depth at this position
        cov_depths.append(cov_depth) # append the new coverage depth to the vector
        if cov_depth == 0: # not covered by any short read
            continue # the base is kept the same as the reference sequence
        else:
            if ref_ins_dict.has_key(i+1): # If there is insertion into the ref sequence at this position, ref_ins_dict is 1-based, unlike the non_ins_bps
                ins_bps = ref_ins_dict[i+1]
                if max(ins_bps) > 0.5 * cov_depth or max(ins_bps) > max(non_ins_bps): # If more than half of the aligned reads have same base pair inserted at this position,\
                                                                                      # or this is the most common vote among all options (substitution and deletion), 
                                                                                      # consider this a 'true' insertion
                    true_ins[i] = alphabet[ins_bps.index(max(ins_bps))]
                    sys.stdout.write('At position {}: inserted {}\t{:>}\n'.format(i + 1, alphabet[ins_bps.index(max(ins_bps))], float(max(ins_bps))/float(cov_depth)))
            else: # no insertion detected, see which one of 'ACGTD' has the highest frequency
                consensus_bp = alphabet[non_ins_bps.index(max(non_ins_bps))]
                if consensus_bp != rSeq[i]: # print message for changed base
                    sys.stdout.write('At position {}: {} --> {} \t{:>}\n'.format(i + 1, rSeq[i], consensus_bp, float(max(non_ins_bps))/float(cov_depth)))
                    rSeq = rSeq[:i] + consensus_bp + rSeq[(i+1):] # change the nucleotide at this position, it could be a 'D' 

    # Now insert the 'true' insertions into the sequence
    length_increase = 0
    for ins in sorted(true_ins):
        insert_position = ins + length_increase
        rSeq = rSeq[:insert_position] + true_ins[ins] + rSeq[insert_position:]
        length_increase += 1 # increased length of the reference sequence

    # Now delete the 'true' deletions from the sequence
    sys.stdout.write('Found {} deletions\n'.format(rSeq.count('D')))
    rSeq = re.sub('D','',rSeq)

    return rSeq, cov_depths

## ======================================================================
## Read a single sequence from fasta file, return the sequence as a string
## ======================================================================
def read_single_seq(fastaFile):
    ''' Read a single sequence from fasta file.
    '''
    with open(fastaFile,'r') as fasta:
        seq = fasta.read().split('\n',1)[-1]
    seq = re.sub('\n','',seq)
    #print seq
    return seq

## ======================================================================
## Find consensus reference sequence from the sam file and the reference sequence fasta file
## ======================================================================
def get_consensus(samFile, ref_fasta, maxSub=3, maxIns=3, maxDel=3, maxErrRate=0.20):
    rSeq = read_single_seq(ref_fasta)
    ref_bps, ref_ins_dict = read_sam(samFile)
    return ref_consensus(ref_bps, ref_ins_dict, rSeq)

## ======================================================================
## compress repeated letters in a sequence (homopolymer), 
## and find the counts for each non-repeat letter
## ======================================================================
def compress_homopolymer(string):
    ''' Compress repeated letters in a string, and find the counts for each of the non-repeat letter. 
        It converts all letters to upper case first.
        Input:  string - a string
        Output: tuple - (nonrepeating letters in upper case, their corresponding counts)
    '''
    if not string.isupper():
        string = string.upper()
    reduced_string = re.sub(r'([ACGT])\1+',r'\1',string) # Find the reduced string
    counts = []
    current_pos = 0
    # find the counts (lengths of stretches) for each letter in the reduced string
    for letter_index in xrange(len(reduced_string)-1):
        next_diff_pos = string.index(reduced_string[letter_index + 1], current_pos)
        counts.append(next_diff_pos - current_pos)
        current_pos = next_diff_pos
    counts.append(len(string) - current_pos)
    return reduced_string, counts

def locate_in_compress(string, pos):
    reduced_string, counts = compress_homopolymer(string)
    if pos == 0:
        return 0
    else:
        for i in xrange(1,len(reduced_string)):
            if sum(counts[:i]) - 1 >= pos:
                return i-1 
## ======================================================================
## Massage the alignement locally, to prefer indels over substitutions,
## without increasing total number of changes of base pairs
## ======================================================================
# First version takes care of left-hand-side part
def massage_mapping(rSeq, qSeq, rstart, rend, cigar_string):
    ''' cigarstring should be sam 1.4 format, with X indicating mismatch
    '''
    if cigar_string.find('X') == -1 : # no mismatch in this read, nothing to do, just return 0
        return 0 

    char = re.findall('\D',cigar_string) # operation characters, MIDNSHP=X
    char = [x.upper() for x in char]  # convert to upper case
    count = map(int,re.findall('\d+',cigar_string))
    cigar_list = re.findall('\D|\d+',cigar_string) # cigar string as a list of str (numbers and operations) 
    alphabet = 'ACGT'

    left_clip = False
    right_clip = False
    chopped_cigar_l  = ''
    chopped_cigar_r  = ''
    # Clipping on the left end, remove the first operation, and remove the clipped part in the qSeq
    if char[0] in 'SH':
        qSeq = qSeq[count[0]:]
        char = char[1:]
        count = count[1:]
        chopped_cigar_l = ''.join(cigar_list[:2])
        left_clip = True

    # Clipping on the right end, remove the first operation, and remove the clipped part in the qSeq
    if char[-1] in 'SH':
        qSeq = qSeq[:count[-1]]
        char = char[:-1]
        count = count[:-1]
        chopped_cigar_r = ''.join(cigar_list[-2:])
        right_clip = True

    # find the index of mismatch operations in the cigar operations
    firstX_in_cigar = ''.join(char).find('X')
    lastX_in_cigar = ''.join(char).rfind('X')
    numX_in_cigar = char.count('X') # total number of Xs in the cigar string

    sys.stdout.write('first mismatch is the {}th operation.\n'.format(firstX_in_cigar))
    sys.stdout.write('last mismatch is the {}th operation.\n'.format(lastX_in_cigar))
    # Only consider the mismatches in the first 2 positions and the last 2 positions
    if firstX_in_cigar not in [0, 1] and lastX_in_cigar not in [len(char)-2, len(char)-1]:
        return 0

    ### Left end operations ###
    elif firstX_in_cigar in [0,1]:
        # First find the position of the first substitution
        # both in the query sequence and the reference sequence
        firstX = sum(count[x] for x in xrange(len(count)) if x < firstX_in_cigar and char[x] != 'D')  # corresponding position on the query sequence for the first mismatch
        firstX_rSeq = sum(count[x] for x in xrange(len(count)) if x < firstX_in_cigar and char[x] != 'I') - 1 + rstart # corresponding position on the reference sequence for the first mismatch
        sys.stdout.write('first mismatch is the {}th position in query sequence.\n'.format(firstX))
        sys.stdout.write('first mismatch is the {}th position in ref sequence.\n'.format(firstX_rSeq))

        # check and see if the mismatched base is the same as the next base in the reference sequence
        # If so, change the mismatch to insert, shift mapping start position 1 bp to the right
        if qSeq[:(firstX+1)] == rSeq[(firstX_rSeq - firstX+1):(firstX_rSeq + 2)]:
            numX_in_cigar -= 1
            cigar_string = re.sub('1X','1I',cigar_string,count=1)
            rstart += 1
            return rstart, cigar_string
        # Next check if the mismatched base is the same as the previous base in the reference sequence
        # If so, change the mismatch to deletion, shift mapping start position 1 bp to the left
        elif rstart - firstX > 1: # there are more bases to shift the query sequence to the left
            if qSeq[:(firstX + 1)] == rSeq[(firstX_rSeq-1-firstX):(firstX_rSeq)]:
                numX_in_cigar -= 1
                cigar_string = chopped_cigar_l + str(count[firstX_in_cigar - 1] + 1) if firstX==0 else '1'  + '=1D' + ''.join(cigar_list[(2*firstX_in_cigar+2):])
                rstart -= 1
                return rstart, cigar_string
            else:
                return 0
        else:
            return 0

    ### Right end operations ###
    elif numX_in_cigar > 0 and lastX_in_cigar in [len(char)-2, len(char)-1]:
        #print count
        #print char
        lastX = sum(count[x] for x in xrange(len(count)) if x < lastX_in_cigar and char[x] != 'D') # corresponding position on the query sequence for the last mismatch
        #print [count[x] for x in xrange(len(count)) if x < lastX_in_cigar and char[x] != 'D']
        lastX_rSeq = sum(count[x] for x in xrange(len(count)) if x < lastX_in_cigar and char[x] != 'I') - 1 + rstart # corresponding position on the reference sequence for the last mismatch
        #print [count[x] for x in xrange(len(count)) if x < lastX_in_cigar and char[x] != 'I']
        #print lastX_in_cigar, lastX, lastX_rSeq
        sys.stdout.write('last mismatch is the {}th position in query sequence.\n'.format(lastX))
        sys.stdout.write('last mismatch is the {}th position in ref sequence.\n'.format(lastX_rSeq))

        right_chars_len = len(qSeq) - lastX # length of the characters to the right of the mismatch position in the query sequence, this position included

        # check and see if the chars from the mismatched position to the end are the same as the stretch of chars starting from the char to the left in the reference sequence
        # If so, change the mismatch to insertion, shift mapping ending position 1 bp to the left, nothing need to do for rstart
        if qSeq[lastX:] == rSeq[(lastX_rSeq - 1):(lastX_rSeq - 1 + right_chars_len)]:
            print 'shift left'
            if left_clip:
                cigar_string = ''.join(cigar_list[:(2*lastX_in_cigar+3)]) + 'I' + ''.join(cigar_list[(2*lastX_in_cigar + 4):])
                return rstart, cigar_string
            else:
                cigar_string = ''.join(cigar_list[:(2*lastX_in_cigar+1)]) + 'I' + ''.join(cigar_list[(2*lastX_in_cigar + 2):])
                return rstart, cigar_string
        # Next check if the mismatched base is the same as the next base in the reference sequence
        # If so, change the mismatch to deletion
        elif lastX_rSeq + right_chars_len < len(rSeq): # there are more bases to shift the query sequence to the right
            if qSeq[lastX:] == rSeq[(lastX_rSeq + 1):(lastX_rSeq + 1 + right_chars_len)]: # can shift to right and change mismatch to deletion
                print 'shift right'
                print cigar_list
                if right_clip: # if the read is right clipped (theoretically there cannot be any extension to the right, since reads are only chopped at the edge of reference sequence, but anyway...)
                    if cigar_list[-3] == 'X': # mismatch is the last operation
                        cigar_string = ''.join(cigar_list[:-4]) + '1D1=' + chopped_cigar_r
                        return rstart, cigar_string
                    elif cigar_list[-3] == '=': # mismatch is the second last operation
                        cigar_string = ''.join(cigar_list[:-6]) + '1D' + str(int(cigar_list[-4]) + 1) + '=' + chopped_cigar_r
                        return rstart, cigar_string
                    else:
                        return 0
                else:
                    if cigar_list[-1] == 'X': # mismatch is the last operation
                        cigar_string = ''.join(cigar_list[:-2]) + '1D1='
                        return rstart, cigar_string
                    elif cigar_list[-1] == '=': # mismatch is the second last operation
                        cigar_string = ''.join(cigar_list[:-4]) + '1D' + str(int(cigar_list[-2]) + 1) + '='
                        return rstart, cigar_string
                    else:
                        return 0
            else: # cannot shift to right and change to deletion
                return 0
        else:
            return 0

def shift_ends(samFile, rSeq, samNew):
    newsam = open(samNew,'w')
    changeRec = 0
    with open(samFile,'r') as mysam:
        for line in mysam:
            if line[0] == '@': # copy header lines
                newsam.write(line)
            else:
                fields = line.strip('\n').split('\t')
                cigar_string = fields[5]
                if cigar_string.find('X') != -1:
                    print cigar_string
                    qSeq = fields[9]
                    rstart = int(fields[3])
                    rend = rstart + cigar(cigar_string)['ref_len'] - 1
                    res = massage_mapping(rSeq, qSeq, rstart, rend, cigar_string)
                    print res, fields[0]
                    if res != 0:
                        fields[3] = str(res[0])
                        fields[5] = res[1]
                        line = '\t'.join(fields) + '\n'
                        sys.stdout.write(fields[0] + '\n')
                        newsam.write(line)
                        changeRec += 1
                    else:
                        newsam.write(line)
#                if fields[0] == 'HISEQ03:379:C2WP8ACXX:7:1101:4141:2993/2':
#                    break
                else:
                    newsam.write(line)
    newsam.close()
    sys.stdout.write('Total number of records changed is {}. \n'.format(changeRec))

## ======================================================================
##  Get the polymorphic positions and the positions where the insertion happened
## ======================================================================
def get_pos(ref_bps, ref_ins_dict):
    ''' Get the polymorphic positions and the positions where the insertion happened
    '''
    # TODO: however, if the coverage depth is very low (for example 2, one should not discard the bases with support of 1 read #
    #num_bases = [ 5 - x.count(0) - x.count(1) for x in ref_bps ] # for each position, get the number of bases with support of at least 2 reads. If there is only 1 read, it most probably will be an error.
    num_bases = [ 5 - x.count(0) for x in ref_bps ] # for each position, get the number of bases with support of at least 2 reads. If there is only 1 read, it most probably will be an error.
    return [ x for x in xrange(len(num_bases)) if num_bases[x] > 1 ], ref_ins_dict.keys()

def ref_extension(ref_bps, ref_ins_dict, rSeq):
    ''' From the output of function get_consensus ( list of mat/mismat/del positions, and dictionary of ins positions), extend the reference sequence to include all the insertion positions
        
        Input: ref_bps - list of lists, one for each position on the reference sequence (without padding)
               ref_ins_dict - dictionary of insertions, one for each insertion position found in the alignments.
               rSeq - original reference sequence

        Output: The extended consensus sequence, the insertion dictionary where insertions were added to the reference sequence, and coverage depth vector
    '''
    alphabet = 'ACGTD' # all the possible base pairs at a position
    rLen = len(rSeq)
    cov_depths = []
    true_ins = dict() # 'true' insertions into reference sequence

    ## correspondence between old positions and new positions after insertion in the reference sequence
    orig_pos_dict = dict()
    ins_pos_dict = dict()

    # loop through all the positions of the reference sequence
    for i in xrange(len(rSeq)):
        non_ins_bps = ref_bps[i] # non-insertion base calls
        cov_depth = sum(non_ins_bps) # coverage depth at this position
        cov_depths.append(cov_depth) # append the new coverage depth to the vector
        if cov_depth == 0: # not covered by any short read
            continue # the base is kept the same as the reference sequence
        else: # if there is coverage, check and see if there is possible insertion at this position
            if ref_ins_dict.has_key(i+1): # If there is insertion into the ref sequence at this position, ref_ins_dict is 1-based, unlike the non_ins_bps
                ins_bps = ref_ins_dict[i+1] # count vector for all possible insertion bases
                if (cov_depth <= 3) or (4 - ins_bps.count(0) - ins_bps.count(1) > 1  and cov_depth > 3): # if coverage is less than or equal to 3, or at least 2 positions with support greater than 1, treat as insertion 
                    true_ins[i] = alphabet[ins_bps.index(max(ins_bps))] # the most frequent base will be saved
                    sys.stdout.write('At position {}: inserted {}\t{:>}\n'.format(i + 1, alphabet[ins_bps.index(max(ins_bps))], float(max(ins_bps))/float(cov_depth)))

    rSeq_ext = rSeq
    # Now insert the 'true' insertions into the sequence
    length_increase = 0
    for ins in sorted(true_ins):
        insert_position = ins + length_increase
        ins_pos_dict[ins] = insert_position
        while pos < insert_position:
            orig_pos_dict[pos] = pos
            pos += 1
        rSeq_ext = rSeq_ext[:insert_position] + true_ins[ins] + rSeq_ext[insert_position:]
        length_increase += 1 # increased length of the reference sequence

    return rSeq_ext, true_ins, orig_pos_dict, ins_pos_dict, cov_depths

def getinfo_at_ambipos(samFile, true_ins, poly_pos, maxSub=3, maxIns=3, maxDel=3, maxSubRate=0.02, maxInsRate=0.2, maxDelRate=0.2):
    ''' Get support information (of short reads) at the polymorphic positions
    '''
    alphabet = 'ACGTD' # all the possible base pairs at a position
    reads_dict = dict()
    proc_reads = 0
    with open(samFile, 'r') as mysam:
        for line in mysam:
            if line[0] == '@': # header line
                continue
            else:
                line = line.strip()
                if not is_record_bad(line,maxSub, maxIns, maxDel, maxSubRate, maxInsRate, maxDelRate): # if this alignment record passes the threshold, gather its information
                    fields = line.split('\t')
                    rstart = int(fields[3])
                    cigar_string = fields[5]
                    rend = rstart + cigar(cigar_string)['ref_len'] - 1
                    pos_dict, ins_dict = get_bases(line) # base calling from this read
                    read_string = ''
                    ins_string = ''

                    for pos in poly_pos: # polymorphic positions
                        if pos < max(pos_dict.keys()):
                            if pos >= rstart-1:
                                read_string += str(pos) + pos_dict[pos+1][-1]
                        else:
                            break

                    for pos in true_ins: # insertion positions
                        if len(ins_dict) == 0:
                            if pos >= rstart - 1:
                                ins_string += str(pos) + 'D'
                        else:
                            if pos < max(ins_dict.keys()):
                                if pos >= rstart - 1 and pos < rend :
                                    if ins_dict.has_key(pos+1):
                                        ins_string += str(pos) + ins_dict[pos+1][-1]
                                    else:
                                        ins_string += str(pos) + 'D'
                            else:
                                break

                    all_string = read_string + ':' + ins_string # read_string and ins_string concatenated
                    if reads_dict.has_key(all_string):
                        reads_dict[all_string] += 1
                    else:
                        reads_dict[all_string] = 1
                    proc_reads += 1

                    if proc_reads % 10000 == 0:
                        sys.stdout.write("processed {} good reads\n".format(proc_reads))
                        return reads_dict
    return reads_dict

## ======================================================================
## use dynamic progamming (Needleman-Wunch) to find all best alignments, 
## the following shift the indel positions to the leftmost position in all the
## output alignment, and pick out the non-equivalent ones
## ======================================================================
def shift_to_left(align):
    ''' Take the result from Bio.pairwise2.align.global**, shift the indels in the homopolymer to the leftmost positions.
    
        Input:  The result (align) is a list of tuples: (seqA, seqB, score, begin, end). seqA and seqB are strings showing the alignment between the sequences. score is the score of the alignment. begin and end are indexes into seqA and seqB that indicate the where the alignment occurs.
                Consider one tuple in the list as input to this function.
                Note: We'll consider that seqA is part of the reference sequence, and seqB is the short read sequence, and the global alignment is done, with no gap penalty at the ends for seqB, but gap penalty at the ends of seqA

        Output: The alignment after shifting
    '''
    seqA, seqB, score, begin, end = align
    # First find all the insertion positions, ignoring the opening and ending gaps in seqB
    first_non_gap = 0 if seqB[0]!='-' else re.search(r'^[-]+',seqB).end()
    last_non_gap = len(seqB) if seqB[-1]!='-' else re.search(r'[-]+$',seqB).start()

    ins_char = re.compile('-')

    # insert positions in sequence A and B
    insA = [m.start() for m in ins_char.finditer(seqA, first_non_gap, last_non_gap)]
    insB = [m.start() for m in ins_char.finditer(seqB, first_non_gap, last_non_gap)]
    insAB = insA + insB
    insAB.sort()

    for ins_pos in insAB:
        if ins_pos in insA:
            base = seqB[ins_pos] # base call corresponding to the insertion position
            if seqA[ ins_pos - 1 ] == base: # if the base to the left is the same, shift it
                l_base = re.search(base+'+$', seqA[:ins_pos]).start()
                seqA = seqA[:l_base] + '-' + seqA[(l_base+1):ins_pos] + base + seqA[(ins_pos+1):]
        else:
            base = seqA[ins_pos] # base call corresponding to the insertion position
            if seqB[ ins_pos - 1 ] == base: # if the base to the left is the same, shift it
                l_base = re.search(base+'+$', seqB[:ins_pos]).start()
                seqB = seqB[:l_base] + '-' + seqB[(l_base+1):ins_pos] + base + seqB[(ins_pos+1):]

    return (seqA, seqB, score, begin, end)

def shift_and_reduce(align_list):
    ''' Take the result from Bio.pairwise2.align.global**, shift the indels in the homopolymer to the leftmost positions, and keep the distinct ones after shifting
        Input:  list of tuples: (seqA, seqB, score, begin, end)
        Output: list of tuples where homopolymer equivalent ones are reduced
    '''
    if len(align_list) == 1:
        return align_list
    else:
        new_list = []
        for align in align_list:
            new_align = shift_to_left(align)
            if new_align not in new_list:
                new_list.append(new_align)
        return new_list
    
def pick_align(align_list):
    ''' From a list of equivalent alignments between 2 sequences, pick the one whose indel positions are the most left
    '''
    leftmost_indel_pos = (10000,10000)
    bestalign = ''
    for align in align_list:
        seqA, seqB, score, begin, end = align
        # First find all the insertion positions, ignoring the opening and ending gaps in seqB
        first_non_gap = 0 if seqB[0]!='-' else re.search(r'^[-]+',seqB).end()
        last_non_gap = len(seqB) if seqB[-1]!='-' else re.search(r'[-]+$',seqB).start()

        # Only look at the aligned part
        seqA = seqA[first_non_gap:last_non_gap]
        seqB = seqB[first_non_gap:last_non_gap]

        ins_char = re.compile('-')
        # insert positions in sequence A and B
        insA = [m.start() for m in ins_char.finditer(seqA)]
        insB = [m.start() for m in ins_char.finditer(seqB)]
        this_indel_pos = (sum(insA+insB), sum(insB)) # use two sums as the measure, one for both sequences and one for the query sequence

        if this_indel_pos < leftmost_indel_pos:
            leftmost_indel_pos = this_indel_pos
            bestalign = align
    return seqA, seqB, score, first_non_gap, last_non_gap # return the list of aligns whose indel positions are the leftmost
    #return bestalign # return the list of aligns whose indel positions are the leftmost

def get_cigar(seqA, seqB):
    ''' Get CIGAR string from the align result 
        Input:  seqA and seqB from the alignment (seqA, seqB, score, begin, end)
        Output: cigar_string, seqA and seqB's aligned part (with opening and ending gaps trimmed)
    '''

    # First find all the insertion positions, ignoring the opening and ending gaps in seqB
    first_non_gap = 0 if seqB[0]!='-' else re.search(r'^[-]+',seqB).end()
    last_non_gap = len(seqB) if seqB[-1]!='-' else re.search(r'[-]+$',seqB).start()
    seqA = seqA[first_non_gap:last_non_gap]
    seqB = seqB[first_non_gap:last_non_gap]
    
    cigar_dict = {'=':0, 'X':0, 'D':0, 'I':0}
    cigar_string = ''
    mode = '0'
    if seqA[0] == seqB[0]:
        mode = '='
        cigar_dict['='] = 1
    elif seqA[0] == '-':
        mode = 'I'
        cigar_dict['I'] = 1
    elif seqB[0] == '-':
        mode = 'D'
        cigar_dict['D'] = 1
    elif seqA[0] != seqB[0]:
        mode = 'X'
        cigar_dict['X'] = 1
    else:
        sys.stderr.write('{} <-> {} in alignment! \n'.format(seqA[0], seqB[0]))

    for i in xrange(1,len(seqA)):
        if seqA[i] == seqB[i]:
            if mode == '=':
                cigar_dict['='] += 1
            else:
                cigar_string += str(cigar_dict[mode]) + mode
                cigar_dict[ '=' ] = 1
                mode = '='
        elif seqA[i] == '-':
            if mode == 'I':
                cigar_dict['I'] += 1
            else:
                cigar_string += str(cigar_dict[mode]) + mode
                cigar_dict['I'] = 1
                mode = 'I'
        elif seqB[i] == '-':
            if mode == 'D':
                cigar_dict['D'] += 1
            else:
                cigar_string += str(cigar_dict[mode]) + mode
                cigar_dict['D'] = 1
                mode = 'D'
        elif seqA[i] != seqB[i]: 
            if mode == 'X':
                cigar_dict['X'] += 1
            else:
                cigar_string += str(cigar_dict[mode]) + mode
                cigar_dict['X'] = 1
                mode = 'X'
    cigar_string += str(cigar_dict[mode]) + mode
    return cigar_string,first_non_gap
## ======================================================================
## Extract the good region of PacBio read, which has >= 10 coverage depth,
## and is 1000 bps long (contiguous).
## If there is more than 1 region in the PacBio read satisfying this condition,
## keep them both, separately (as in 2 files??)
## ======================================================================
def get_good_regions(ref_bps, ref_ins_dict, rSeq, minPacBioLen=1000, minCV=10):
    ''' From mapping results from Illumina reads to a PacBio read, find the regions of PacBio read that satisfy these conditions:
        1.  Every base pair in the region is covered by at least 10 reads
        2.  The contiguous region has length >= 1000 bps
        
        Input:  ref_bps - dictionary for the match/mismatch and the deletion positions
                ref_ins_dict - dictionary for the insertion positions
                rSeq - PacBio sequence
        Output: ? TODO: fill in this part
    '''
    alphabet = 'ACGTD' # all the possible base pairs at a position
    rLen = len(rSeq)
    # if the whole sequence is shorter than 1000 bps, do not consider this sequence
    if rLen < minPacBioLen:
        return []

    cov_depths = []
    good_regions = []
    # loop through all the positions of the reference sequence
    for i in xrange(rLen):
        non_ins_bps = ref_bps[i] # non-insertion base calls
        cov_depth = sum(non_ins_bps) # coverage depth at this position
        cov_depths.append(cov_depth) # append the new coverage depth to the vector
    
    low_CV_pos = [-1] + [ i for i in xrange(rLen) if cov_depths[i] < minCV ] + [rLen]
    for i in xrange(1, len(low_CV_pos)):
        if low_CV_pos[i] - low_CV_pos[i-1] >= minPacBioLen:
            good_regions.append((low_CV_pos[i-1]+1, low_CV_pos[i])) # found a good region (begin, end) where rSeq[begin:end] is long enough and covered well
    return good_regions
