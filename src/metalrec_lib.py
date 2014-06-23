#!/usr/bin/python

''' Collection of functions used in MetaLREC project
'''
import re
import sys

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
## parse a single line in sam file with an aligned read record
## extract information
## This is somehow redundant with the samread class implementation.
## Question: Will this be faster than construct a whole class??
## ======================================================================
def readAlign(alignRecord):
    """ From alignment record for a read, get the following information:
        1. qname - query read name
        2. flag
        3. rname - reference sequence name
        4. pos - 1-based leftmost mapping position on the reference sequence
        5. mapQ - mapping quality
        6. cigarstring - CIGAR string
        the files following are: rnext, pnext, tlen, seq, qual
        10. seq - segment sequence
        
        Output:
        Dictionary with these fields, including some other inferred ones:
        1. qstart, qend, rstart(pos), rend, positions
        2. what else?
    """
    fields = alignRecord.split("\t") # split by tabs
    qname = fields[0]
    flag = int(fields[1]) # bitwise Flag

    # FLAG explanation one by one
    is_paired = flag & 0x1 == 0x1 # 0x1: read is paired
    is_proper_pair = flag & 0x2 == 0x2 # 0x2: reads mapped in proper pair 
    is_unmapped = flag & 0x4 == 0x4 # 0x4: read unmapped
    mate_is_unmapped = flag & 0x8 == 0x8 # 0x8: mate/next read unmapped
    is_reverse = flag & 0x10 == 0x10 # 0x10: read being reverse complemented
    mate_is_reverse = flag & 0x20 == 0x20 # 0x20: mate/next read being reverse complemented
    is_read1 = flag & 0x40 == 0x40 # 0x40: first read in the pair
    is_read2 = flag & 0x80 == 0x80 # 0x80: second read in the pair
    is_secondary = flag & 0x100 == 0x100 # 0x100: secondary alignment, false if mapping is primary
    is_qcfail = flag & 0x200 == 0x200 # 0x200: not passing quality control
    is_duplicate = flag & 0x400 == 0x400 # 0x400: read is PCR or optical duplicate
    is_supplementary = flag & 0x800 == 0x800 # 0x800: supplementary alignment (part of a chimeric alignment)
    
    rname = fields[2] # reference sequence name
    rstart = int(fields[3]) # starting mapping position on the reference sequence, 1-based
    try: # mapping quality
        mapQ = int(fields[4])
    except ValueError:
        mapQ = -1
    cigarstring = fields[5] # CIGAR string
    qSeq = fields[9] # segment sequence

    return {'qname':qname, 'flag':flag,'is_paired':is_paired, 'is_proper_pair':is_proper_pair, 'is_unmapped':is_unmapped, 'mate_is_unmapped': mate_is_unmapped, 'is_reverse':is_reverse, 'mate_is_reverse':mate_is_reverse, 'is_read1':is_read1, 'is_read2':is_read2, 'is_secondary':is_secondary, 'is_qcfail':is_qcfail, 'is_duplicate':is_duplicate, 'is_supplementary':is_supplementary, 'rname':rname, 'rstart':rstart, 'mapQ':mapQ, 'cigarstring':cigarstring, 'qSeq':qSeq }

## ======================================================================
## function to check the CIGAR string to see if this record is bad or not
## works with sam format 1.4
## later should make it work with sam format 1.3
## ======================================================================
def is_record_bad(alignRecord,maxSub=3, maxIns=3, maxDel=3,maxErrRate=0.20):
    ''' Test and see if an alignment record is bad in the sam file.
        Input:  alignRecord - a mapping line in sam file
                maxSub - maximum stretches of substitution
                maxIns - maximum stretches of insertion
                maxDel - maximum stretches of deletion
                maxErrRate - maximum total error rate of the read
        Output: boolean value, True if bad else False
    '''
    fields = alignRecord.split("\t") # split by tabs
    cigarstring = fields[5] # CIGAR string
    if cigarstring == '*': # if cigar string is not available, treat as good unless the tag indicates that the read is unmapped
        if readAlign(alignRecord)['is_unmapped']:
            return True
        else:
            return False
    else:
        cigar_info = cigar(cigarstring)
        # if consecutive sub or indels is longer than the threshold, treat as bad
        if cigar_info['max_ins'] > maxIns or cigar_info['max_sub'] > maxSub or cigar_info['max_del'] > maxDel:
            return True
        else:
            mismatchLen = cigar_info['ins_len'] + cigar_info['del_len'] + cigar_info['sub_len']
            # if total error bps is too big, treat as bad
            if mismatchLen > (cigar_info['seq_len'] + cigar_info['del_len'])*maxErrRate:
                return True
            else:
                return False
    
## ======================================================================
## Remove low quality alignments from short reads to long read, from the 
## mapping results.
## ======================================================================
def rm_bad_record(samFile,samNew, maxSub=3, maxIns=3, maxDel=3, maxErrRate=0.20):
    newsam = open(samNew,'w')
    keepRec = 0
    with open(samFile,'r') as mysam:
        for line in mysam:
            if line[0] == '@': # copy header lines
                newsam.write(line)
            else:
                record = line.strip('\n')
                if not is_record_bad(record, maxSub, maxIns, maxDel, maxErrRate): # if this alignment is good
                    newsam.write(line)
                    keepRec += 1
    newsam.close()
    sys.stdout.write('Total number of records kept is {}. \n'.format(keepRec))

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

## ======================================================================
## From the mapping sam file and the reference sequence,
## get the mapping information for the reference sequence, base by base
## ======================================================================
def read_sam(samFile, maxSub=3, maxIns=3, maxDel=3, maxErrRate=0.20):
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
                if not is_record_bad(line,maxSub, maxIns, maxDel, maxErrRate): # if this alignment record passes the threshold, gather its information
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

def getinfo_at_ambipos(samFile, true_ins, poly_pos, maxSub=3, maxIns=3, maxDel=3, maxErrRate=0.20):
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
                if not is_record_bad(line,maxSub, maxIns, maxDel, maxErrRate): # if this alignment record passes the threshold, gather its information
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
