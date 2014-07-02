#!/usr/bin/python

''' Collection of functions used in MetaLREC project
'''
import re
import sys
import samread
from Bio import pairwise2 # pairwise alignment using dynamic programming
from Bio.pairwise2 import format_alignment

## ======================================================================
def dict_to_string(dictionary):
    ''' Convert dictionary to a string, key value all concatenated. Dictionary is first sorted by keys '''
    mylist = []
    for key in sorted(dictionary):
        mylist.append(str(key))
        mylist.append(str(dictionary[key]))
    return ''.join(mylist)
## ======================================================================
def cigar(cigar_string):
    """ parse CIGAR string from .sam file, find number of matches, mismatches, insertion (to ref), deletion,
        soft clipping, hard clipping, padding, maximum continuous substitution/insertion/deletion length
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
def md(MD_tag):
    """
    Given MD tag and a sequence, find the number of matched, deletion, and substitution bps 

    Input:  MD tag as a string
    Output: 1-level dictionary with matching, deletion, and substitution lengths
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
def read_single_seq(fastaFile):
    ''' Read a single sequence from fasta file.
        Input:  fastaFile - fasta file with 1 sequence in it
        Output: seq - sequence as a string
    '''
    with open(fastaFile,'r') as fasta:
        seq = fasta.read().split('\n',1)[-1]
    seq = re.sub('\n','',seq)
    #print seq
    return seq
## ======================================================================
def is_record_bad(alignRecord,maxSub=3, maxIns=3, maxDel=3,maxSubRate=0.02, maxInsRate=0.2, maxDelRate=0.2):
    ''' Test and see if an alignment record is bad in the sam file. 
        TODO: only works with sam v1.4 cigar string now, combine the ref seq to also work for v1.3 string later.
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
        elif cigar_info['sub_len'] > maxSubRate * cigar_info['seq_len'] or cigar_info['ins_len'] > maxInsRate * cigar_info['seq_len'] or cigar_info['del_len'] > maxDelRate * cigar_info['seq_len']:
            return True
        # Finally, if it passes all the thresholds, it's a good record
        else:
            return False
## ======================================================================
def clean_samfile(samFile,samNew, rseq, maxSub=3, maxIns=3, maxDel=3,maxSubRate=0.02, maxInsRate=0.2, maxDelRate=0.2):
    ''' Remove low quality alignments from short reads to long read, from the 
        mapping results, and write reads with good mapping quality to a new sam file for visualization.
        Input:  samFile - sam file from mapping multiple short reads to a long read
                samNew - file name for the reduced good quality mappings
                rseq - reference sequence as a string, used to remap the reads with new penalty
                maxSub, maxIns, maxDel, maxSubRate, maxInsRate, maxDelRate are used in is_record_bad
        Output: samNew - new sam file for visualization of cleaner mapping
    '''
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
                    ref_region_start = max( myread.rstart - 5, 1)
                    ref_region_end = min(myread.get_rend() + 5, rLen)
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
def get_bases(cigar_string, qseq='', start_pos=''):
    ''' from CIGAR string, query segment (in alignment record), and starting position (1-based) on the ref sequence, return position wise base call from the read.
        Assuming that cigar_string is available, doesn't matter if it's sam 1.3 format or 1.4 format
        If there is only one input argument, treat it as the whole alignment record.
        
        Input:  cigar_string - cigar string 
                qseq - query aligned segment (optional)
                start_pos - 1-based starting mapping position on the reference sequence (optional)
        
        Output: (pos_dict, ins_dict), tuple of 2 dictionaries, one for the non-insertion positions, and one for the insertion positions.
                When there is a deletion from the reference sequence, the base called will be "D". The bases are 0-based.
                pos_dict: ref_pos => (query_pos, base), or ref_pos => 'D'
                ins_dict: ref_pos => (query_pos, base(s)) # the inserted length could be 1, or greater than 1
    '''
    # the whole alignment record
    if qseq == '' and start_pos == '':
        fields = cigar_string.strip('\n').split('\t')
        cigar_string = fields[5]
        qseq = fields[9]
        start_pos = int(fields[3])
    # cigar operations and the corresponding counts
    char = re.findall('\D',cigar_string) # operation characters, MIDNSHP=X
    char = [x.upper() for x in char]  # convert to upper case
    count = map(int,re.findall('\d+',cigar_string)) # corresponding count for each operation

    # output dictionaries
    pos_dict = dict() # ref_pos (1-based) => base call
    ins_dict = dict() # ref_pos (1-based) => inserted base

    query_pos = 0 # tracking the position currently under investigation 
    start_pos -= 1 # convert starting position to 0-based instead of 1-based as used in sam file

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
                pos_dict[start_pos] = qseq[query_pos] # base at the ref position in the read
                start_pos += 1
                query_pos += 1
        elif char[pos] == 'I': # insertion into reference sequence
            ins_dict[start_pos] = qseq[query_pos:(query_pos + count[pos])]
            query_pos += count[pos]
        elif char[pos] == 'D' or char[pos] == 'N': # deletion or skipped region from reference sequence
            for i in xrange(count[pos]):
                pos_dict[start_pos] = 'D'
                start_pos += 1
        else: # unknown CIGAR operations, print message
            sys.stdout.write(" unknown CIGAR operation: {}\n".format(char[pos]))

    return pos_dict, ins_dict
## ======================================================================
def get_bases_from_align(align, start_pos):
    ''' Same as get_bases, but from pairwise alignment results (of Bio.pairwise2 module) instead of a cigar string

        Input:  align - alignment output from Bio.pairwise2.globalXX, tuple (seqA, seqB, score, begin, end)
                start_pos - starting position (1-based from sam file) on the reference sequence

        Output: (pos_dict, ins_dict), tuple of 2 dictionaries, one for the non-insertion positions, and one for the insertion positions.
                When there is a deletion from the reference sequence, the base called will be "D". The bases are 0-based.
                pos_dict: ref_pos => (query_pos, base), or ref_pos => 'D'
                ins_dict: ref_pos => (query_pos, base(s)) # the inserted length could be 1, or greater than 1
    '''
    pos_dict = dict() # ref_pos (0-based) => base call
    ins_dict = dict() # ref_pos (0-based) => inserted base
    
    seqA, seqB, score, begin, end = align
    ref_pos = start_pos - 1 # convert 1-based position to 0-based position

    for i in xrange(len(seqA)):
        if seqA[i] == '-': # insertion
            if ref_pos in ins_dict: # insertion of more than 1 base pair
                ins_dict[ ref_pos ] = ins_dict[ref_pos] + seqB[i]
            else:
                ins_dict[ ref_pos ] = seqB[i]
        elif seqB[i] == '-': # deletion
            pos_dict[ ref_pos ] = 'D'
            ref_pos += 1
        else: # match or mismatch
            pos_dict[ ref_pos ] = seqB[i]
            ref_pos += 1

    return pos_dict, ins_dict
## ======================================================================
def read_and_process_sam(samFile,rseq, maxSub=3, maxIns=3, maxDel=3,maxSubRate=0.02, maxInsRate=0.2, maxDelRate=0.2, minPacBioLen=1000, minCV=10):
    ''' Get consensus sequence from alignments of short reads to a long read, in the process, filter out bad reads and improve mapping

        Input:  samFile - sam file generated by mapping
                rseq - reference sequence as a string
                maxSub, maxIns, maxDel, maxSubRate, maxInsRate, maxDelRate are used to filter out badly mapped reads
                minPacBioLen - contiguous region length threshold to be a good region
                minCV - minimum coverage depth for a position to be considered as part of a good region

        Output: ref_bps - list of lists, one for each position on the reference sequence (without padding)
                ref_ins_dict - dictionary of insertions, one for each insertion position found in the alignments.
    '''
    alphabet = 'ACGTD' # all the possible base pairs at a position
    keepRec = 0
    discardRec = 0
    lineNum = 0
    
    with open(samFile, 'r') as mysam:
        for line in mysam:
            lineNum += 1 
            if line[0] == '@': # header line
                if line[1:3] == 'SQ': # reference sequence dictionary
                    rname = line[(line.find('SN:') + len('SN:')) : line.find('\t',line.find('SN:'))] # reference sequence name
                    rLen = int(line[(line.find('LN:') + len('LN:')) : line.find('\t',line.find('LN:'))]) # reference sequence length
                    sys.stdout.write("reference sequence {} has length {}. \n".format(rname, rLen))
                    if rLen < minPacBioLen:
                        sys.stdout.write("reference length is smaller than threshold {}.\n".format(minPacBioLen))
                        return 0 
                    ref_bps = [ [0] * 5 for x in xrange(rLen) ]  # list of lists, one list corresponding to each position on the reference sequence
                    ref_ins_dict = dict() # global insertion dictionary for the reference sequence
                    readinfo = dict() # dictionary storing read information (base call for the mapped read)
            else:
                record = line.strip('\n')
                fields = record.split('\t')
                if not is_record_bad(record, maxSub, maxIns, maxDel, maxSubRate, maxInsRate, maxDelRate): # if this alignment is good
                    keepRec += 1
                    myread = samread.SamRead(record)
                    # improve the alignment to the reference sequence by dynamic programming
                    # extend original starting and ending positions on the reference sequence of the mapped region by 5 bps on each end
                    # TODO: change length of extension to consider to a variable, instead of a fixed number 5
                    ref_region_start = max( myread.rstart - 5, 1)
                    ref_region_end = min(myread.get_rend() + 5, rLen)
                    # query sequence with clipped part trimmed
                    trimmed_qseq = myread.get_trim_qseq()
                    # redo global alignment using dynamic programming, indel penalty -0.9, mismatch penalty -1, opening and ending gaps at the query sequence has no penalty
                    realign_res = pairwise2.align.globalms(rseq[(ref_region_start-1):ref_region_end], trimmed_qseq, 0, -1, -0.9, -0.9, penalize_end_gaps=[True, False])
                    new_align = pick_align(realign_res) # pick the best mapping: indel positions are the leftmost collectively
                    pos_dict, ins_dict = get_bases_from_align(new_align, ref_region_start + new_align[3])
                    read_string = dict_to_string(pos_dict) + ':' +  dict_to_string(ins_dict)
                    if read_string in readinfo:
                        readinfo[read_string] += 1
                    else:
                        readinfo[read_string] = 1

                    for pos in pos_dict: # all the matching/mismatching/deletion positions
                        ref_bps[pos][alphabet.find(pos_dict[pos])] += 1 # update the corresponding base call frequencies at the position

                    # Insertion positions dictionary. Every position where insertion happens has a dictionary, with position mapped to list of length 4 lists
                    # In this list, every base pair has its list, if only 1 char is inserted, then the list has length 1
                    # For example, {3:[[0,10,0,9]]} means there is only 1 bp insertion for all the reads, 10 inserted 'C' and 9 inserted 'T'
                    #              {3:[[0,10,0,9], [1,10,0,0]] means there are 2bps insertion, for the second char, 1 read inserted 'A', 10 reads inserted 'C'
                    for ins in ins_dict: # all the insertion positions
                        ins_chars = ins_dict[ins]
                        if ins not in ref_ins_dict: # if this position has not appeared in the big insertion dictionary, initialize it
                            ref_ins_dict[ins] = [[0,0,0,0] for ii in xrange(len(ins_chars))] # length is 4 because there is no 'D' at insertion position
                        elif len(ins_chars) > len(ref_ins_dict[ins]): # if the inserted bps is longer than the list corresponding to this position, make it longer
                            ref_ins_dict[ins] += [[0,0,0,0] for ii in xrange(len(ref_ins_dict[ins]),len(ins_chars))]
                        for i in xrange(len(ins_chars)):
                            #print ins_chars, ref_ins_dict[ins]
                            ref_ins_dict[ins][i][alphabet.find(ins_chars[i])] += 1
                            #print ref_ins_dict[ins]

                    if keepRec % 1000 == 0:
                        sys.stdout.write('  processed {} good records\n'.format(keepRec))
                    # TODO: use only 5000 reads, for testing purpose, remove this for real analysis
                    if keepRec % 5000 == 0:
                        print lineNum
                        sys.stdout.write("discarded {} reads so far.\n".format(discardRec))
                        return ref_bps, ref_ins_dict
                else:
                    discardRec += 1
    
    sys.stdout.write("discarded {} reads.\n".format(discardRec))
    return ref_bps, ref_ins_dict, readinfo

## ======================================================================
## use dynamic progamming (Needleman-Wunch) to find all best alignments, 
## the following shift the indel positions to the leftmost position in all the
## output alignment, and pick out the non-equivalent ones
## ======================================================================
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

## ======================================================================
##  Get the polymorphic positions in the specified region (>=1000bps and covered by >=10 reads)
## ======================================================================
def get_poly_pos(ref_bps, ref_ins_dict, region=None, minReads=3, minPercent=0.01):
    ''' Get the polymorphic positions and the positions where the insertion happened
    '''
    if region is None: # region to consider
        region = (0, len(ref_bps))
    poly_bps = [] # list of tuples
    poly_ins = []
    consensus_bps = []
    consensus_ins = []
    alphabet = 'ACGTD' # all the possible base pairs at a position
    cvs = []
    for pos in xrange(region[0], region[1]):
        # check the match/mismatch/deletion first
        cv = sum(ref_bps[pos]) # coverage depth
        cvs.append(cv)
        base_calls = [ alphabet[i] for i in xrange(5) if ref_bps[pos][i] >= minReads and ref_bps[pos][i] >= cv * minPercent ]
        if len(base_calls) > 1:
            poly_bps.append((pos, base_calls))
        elif len(base_calls) == 1:
            consensus_bps.append((pos, base_calls[0]))
        # check the insertion
        if pos in ref_ins_dict: 
            base_calls = [ alphabet[i] for i in xrange(4) if ref_ins_dict[pos][i] >= minReads and ref_ins_dict[pos][i] >= cv * minPercent ]
            if len(base_calls) > 1:
                poly_ins.append((pos, base_calls))
            elif len(base_calls) == 1:
                consensus_ins.append((pos, base_calls[0]))

    return poly_bps, poly_ins, consensus_bps, consensus_ins


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

