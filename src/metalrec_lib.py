#!/usr/bin/python

''' Collection of functions used in MetaLREC project
'''
import re
import sys, os
sys.path.append("/lustre/atlas/scratch/chaij1/csc124/biopython-1.64/lib/python2.7/site-packages") # for loading modules on titan
from numpy import *
from Bio import pairwise2 # pairwise alignment using dynamic programming
from Bio.pairwise2 import format_alignment
# import local module
import mycolor # print with color in terminal
import samread # for manipulating sam record
alphabet = 'ACGTD'
def minCover(cv):
    ''' Get minimum read support for a base call to be considered correct
        Input:  cv - coverage depth at a position
        Output: min_cover - minimum read support
    '''
    if cv <= 3:
        return 1
    else:
        return 2
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
    pad_len = sum(count[i] for i in xrange(len(char)) if char[i] == 'P') # padded bps, inserted both in reference sequence and the read sequence)
    left_clip_len = 0
    right_clip_len = 0
    if char[0] in 'SH':
        left_clip_len = count[0]
    if char[-1] in 'SH':
        right_clip_len = count[-1]
    clip_len = left_clip_len + right_clip_len # clipped bps, including hard and soft clipping, on both end of the reads

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
        align_len = match_len + sub_len # aligned bps(including matching and substitution)
    else:   # sam 1.3, aligned bps from 'M'
        align_len = sum(count[i] for i in xrange(len(char)) if char[i] == 'M') # aligned bps(including matching and substitution)

    return {'read_len':read_len, 'seq_len': seq_len, 'ref_len':ref_region_len, 'ins_len':ins_len, 'del_len':del_len, 'match_len':match_len, 'sub_len':sub_len, 'align_len':align_len, 'pad_len':pad_len, 'max_ins': max_ins_len, 'max_del': max_del_len, 'max_sub': max_sub_len, 'left_clip_len':left_clip_len, 'right_clip_len':right_clip_len}
## ======================================================================
def md(MD_tag):
    """
    Given MD tag and a sequence, find the number of matched, deletion, and substitution bps 

    Input:  MD tag as a string
    Output: 1-level dictionary with matching, deletion, and substitution lengths
    """
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
def read_fasta(fasta_file,trim=True,reverse=False):
    ''' Read .fasta file, and return dictionary object: seqID => Seq
        Input:
        1. fasta_file: input fasta file
        2. trim: trim ID at first space
        3. reverse: if true, seq => ID, instead of ID => seq
        Output
        1. dictionary object
    '''
    fasta = open(fasta_file,'r')
    ## all the contigs in a list
    contigs = fasta.read().split('\n>')
    fasta.close()
    
    contigs[0] = contigs[0][1:] # trim '>' for the first contig
    # input contigs info: contigID => sequence for this contig
    input_info = dict()
    for contig in contigs:
        if trim:
            contigID = contig.split()[0] # ID of this contig
        else:
            contigID = contig.split("\n",1)[0] # full header
        contigSeq = contig.split("\n",1)[1] # its corresponding sequence
        contigSeq = re.sub("\r","",contigSeq) # remove windows style EOF
        contigSeq = re.sub("\n","",contigSeq) # remove unix style EOF

        if reverse: # contigSeq => contigID
            # if two proteins have the same sequences, store the corresponding IDs in a list
            if input_info.has_key(contigSeq):
                input_info[contigSeq].append(contigID)
            else:
                input_info[contigSeq] = [contigID]
        else:
            input_info[contigID] = contigSeq # link them together in the dictionary
    return input_info
## ======================================================================
def read_single_seq(fastaFile):
    ''' Read a single sequence from fasta file.
        Input:  fastaFile - fasta file with 1 sequence in it
        Output: seq - sequence as a string
    '''
    with open(fastaFile,'r') as fasta:
        seq = fasta.read().split('\n',1)[-1]
    seq = re.sub('\n','',seq) # remove both unix and windows style newline characters
    seq = re.sub('\r','',seq)
    return seq
## ======================================================================
def is_record_bad(alignRecord, maxSub=-1, maxIns=-1, maxDel=-1,maxSubRate=0.05, maxInDelRate=0.3):
    ''' Test and see if an alignment record is bad in the sam file. 
        TODO: only works with sam v1.4 cigar string now, combine the ref seq to also work for v1.3 string later.
        Input:  alignRecord - a mapping line in sam file
                maxSub - maximum stretches of substitution
                maxIns - maximum stretches of insertion
                maxDel - maximum stretches of deletion
                maxSubRate - maximum total substitution rate of the read
                maxInDelRate - maximum total indel rate of the read (put InsRate and DelRate restrictions together into this one)
        Output: boolean value, True if bad else False
    '''
    fields = alignRecord.split("\t") # split by tabs
    if len(fields) < 11: # sam record has at least 11 mandatory fields
        sys.stderr.write('alignRecord is bad: \n {} \n'.format(alignRecord))
        return True
    try: # flag
        flag = int(fields[1])
    except ValueError:
        sys.stderr.write('alignRecord is bad: \n {} \n'.format(alignRecord))
        return True
    try: # cigar string
        cigarstring = fields[5] # CIGAR string
    except IndexError:
        sys.stderr.write('alignRecord is bad: \n {} \n'.format(alignRecord))
        return True

    if flag & 0x4 == 0x4: # read is unmapped (most reliable place to tell if a read is unmapped)
        return True
    elif cigarstring == '*': # if cigar string is not available, treat as bad 
        return True
    else: # read is NOT unmapped and cigar string is available
        cigar_info = cigar(cigarstring)
        indel_len = cigar_info['ins_len'] + cigar_info['del_len']
        # if consecutive sub or indels is longer than the threshold, treat as bad
        if (maxIns != -1 and cigar_info['max_ins'] > maxIns) or (maxSub != -1 and cigar_info['max_sub'] > maxSub) or (maxDel != -1 and cigar_info['max_del'] > maxDel):
            return True
        # substitution rate is relative to the aligned region length on the PacBio sequence
        # indelRate is relative to the total length of the mapped segment of the Illumina sequence
        elif cigar_info['sub_len'] > maxSubRate * cigar_info['seq_len'] or indel_len > maxInDelRate * cigar_info['seq_len']: # subRate * mapped ref region; indelRate * seqLen 
            #print "sub_len: ", cigar_info['sub_len'], ", indel len: ", indel_len #DEBUG
            return True
        # Finally, if it passes all the thresholds, it's a good record
        else:
            return False
## ======================================================================
def get_bases(cigar_string, qseq='', start_pos=None):
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
    if qseq == '' or start_pos is None:
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

        Output: (pos_dict, ins_dict), tuple of 2 dictionaries, one for the non-insertion positions, and one for the insertion positions. 0-based index
                When there is a deletion from the reference sequence, the base called will be "D". The bases are 0-based.
                pos_dict: ref_pos => (query_pos, base), or ref_pos => 'D'
                ins_dict: ref_pos => (query_pos, base(s)) # the inserted length could be 1, or greater than 1
    '''
    pos_dict = dict() # ref_pos (0-based) => base call
    ins_dict = dict() # ref_pos (0-based) => inserted base
    
    seqA, seqB, score, begin, end = align
    if len(seqA) != len(seqB): # sanity check
        sys.stderr.write("get_bases_from_align: seqA and seqB do not have the same length\n{}\n{}".format(seqA, seqB))
    ref_pos = start_pos - 1 # convert 1-based position to 0-based position

    for i in xrange(len(seqB)): # seqA: reference, seqB: read
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
def read_and_process_sam_samread(samFile,rseq, maxSub=-1, maxIns=-1, maxDel=-1,maxSubRate=0.05, maxInDelRate=0.3, minPacBioLen=1000, checkEnds = True, outDir=None, verbose=False):
    ''' Get consensus sequence from alignments of short reads to a long read, in the process, filter out bad reads and improve mapping
        Uses SamRead class instead of calling all the functions

        Input:  samFile - sam file generated by mapping
                rseq - reference sequence as a string
                maxSub, maxIns, maxDel, maxSubRate, maxInsRate, maxDelRate are used to filter out badly mapped reads
                minPacBioLen - contiguous region length threshold to be a good region
                outDir - output directory where the fasta file including good Illumina reads and cleaner sam file will be stored, only in verbose mode
                verbose - switch of verbosity

        Output: ref_bps - list of lists, one for each position on the reference sequence (without padding)
                ref_ins_dict - dictionary of insertions, one for each insertion position found in the alignments.
                readinfo - dictionary of read information, info_string => list of names of reads whose mapped segment corresponds to this info_string
    '''
    alphabet = 'ACGTD' # all the possible base pairs at a position
    keepRec = 0
    discardRec = 0
    lineNum = 0
    rname = ''
    # print verbose information, including reduced alignment and reads that are aligned well
    # make sure the file names are set
    if verbose:
        if outDir is None:
            outDir = os.path.dirname(os.path.abspath(samFile)) + '/EC/'
        else:
            outDir = os.path.abspath(outDir) + '/'
        if not os.path.exists(outDir):
            os.makedirs(outDir)

        outsamFile = outDir + 'scrub.sam'
        sys.stdout.write("write new alignment in sam file {}.\n".format(outsamFile))
        
        outFastaFile = outDir + 'goodreads.fasta'
        sys.stdout.write("write aligned sequences in fasta file {}.\n".format(outFastaFile))

        newsam = open(outsamFile,'w')
        outFasta = open(outFastaFile, 'w')

    ref_bps = []
    ref_ins_dict = dict() # global insertion dictionary for the reference sequence
    readinfo = dict() # dictionary storing read information (base call for the mapped read)
    rLen = len(rseq)
    ref_bps = [ [0] * 5 for x in xrange(rLen) ]  # list of lists, one list(length 5) corresponding to each position on the reference sequence
    header_written = False
    with open(samFile, 'r') as mysam:
        for line in mysam:
            lineNum += 1 
            #####################################
            # header line
            if line[0] == '@': 
                if verbose:
                    newsam.write(line) # if new reduced sam file is required, write the head lines into the new sam file
                if line[1:3] == 'SQ': # reference sequence dictionary
                    rname = line[(line.find('SN:') + len('SN:')) : line.find('\t',line.find('SN:'))] # reference sequence name
                    rLen = int(line[(line.find('LN:') + len('LN:')) : line.find('\t',line.find('LN:'))]) # reference sequence length
                    sys.stdout.write("Sequence name: {} \nLength: {}\n".format(rname, rLen))
                    if rLen < minPacBioLen:
                        if verbose:
                            sys.stdout.write("Reference is shorter than threshold {}.\n".format(minPacBioLen))
                        return 0 
            #####################################
            # mapping record lines
            else:
                record = line
                if ' ' in line:
                    myread = samread.BlasrRead(record)
                else:
                    myread = samread.SamRead(record)
                    #print myread.qname # DEBUG
                if rname == '':
                    rname = myread.rName
                if verbose and not header_written: # write header lines for the scrubbed sam file
                    newsam.write("@HD\tVN:1.4\tSO:unsorted\n")
                    newsam.write("@SQ\tSN:{}\tLN:{}\n".format(rname, rLen))
                    newsam.write("@RG\tID:1\n")
                    newsam.write("@PG\tID:metalrec\n")
                    header_written = True
                if not myread.is_read_bad(maxSub, maxIns, maxDel, maxSubRate, maxInDelRate): # if this alignment is good
                    #sys.stdout.write("realign\n") # DEBUG
                    pos_dict, ins_dict = myread.re_align(maxSub, maxIns, maxDel, maxSubRate, maxInDelRate) # realign read to PacBio sequence
                    if len(pos_dict) + len(ins_dict) > 0:
                        keepRec += 1
                        if verbose:
                            newsam.write(myread.generate_sam_record())
                            outFasta.write('>{}\n{}\n'.format(myread.qname, re.sub('-', '', myread.qSeq)))

                        # update string dictionary for the read information
                        # TODO: currently, exactly same reads are processed multiple times, computation can be reduced
                        read_string = dict_to_string(pos_dict) + ':' +  dict_to_string(ins_dict)
                        # single end reads, different reads don't have same names
                        if read_string not in readinfo: # new read_string(key) in the dictionary
                            readinfo[read_string] = []
                        readinfo[read_string].append(myread.qname)

                        for pos in pos_dict: # all the matching/mismatching/deletion positions
                            ref_bps[pos][alphabet.find(pos_dict[pos])] += 1 # update the corresponding base call frequencies at the position

                        for ins in ins_dict: # all the insertion positions
                            ins_chars = ins_dict[ins]
                            if ins not in ref_ins_dict: # if this position has not appeared in the big insertion dictionary, initialize it
                                ref_ins_dict[ins] = [[0,0,0,0] for ii in xrange(len(ins_chars))] # length is 4 because there is no 'D' at insertion position
                            elif len(ins_chars) > len(ref_ins_dict[ins]): # if the inserted bps is longer than the list corresponding to this position, make it longer
                                ref_ins_dict[ins] += [[0,0,0,0] for ii in xrange(len(ref_ins_dict[ins]),len(ins_chars))]
                            for i in xrange(len(ins_chars)):
                                ref_ins_dict[ins][i][alphabet.find(ins_chars[i])] += 1

                        if verbose and keepRec % 1000 == 0:
                            sys.stdout.write('  processed {} good records\n'.format(keepRec))
                    else:
                        discardRec += 1
                else:
                    discardRec += 1
    
    if verbose:
        newsam.close()
        outFasta.close()
        sys.stdout.write("discarded {} reads, kept {} reads.\n".format(discardRec, keepRec))
        sys.stdout.write("number of unique reads: {}\n".format(len(readinfo)))
    
    return ref_bps, ref_ins_dict, readinfo
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
            look_pos = ins_pos
            while seqA[ look_pos - 1] == '-': # look for the base call to the left, ignore the insertions since they don't matter
                look_pos -= 1
            if seqA[ look_pos - 1 ] == base: # if the base to the left is the same, shift it
                l_base = re.search(base+'+$', seqA[:look_pos]).start()
                seqA = seqA[:l_base] + '-' + seqA[(l_base+1):ins_pos] + base + seqA[(ins_pos+1):]
        else:
            base = seqA[ins_pos] # base call corresponding to the insertion position
            look_pos = ins_pos
            while seqB[ look_pos - 1] == '-':
                look_pos -= 1
            if seqB[ look_pos - 1 ] == base: # if the base to the left is the same, shift it
                l_base = re.search(base+'+$', seqB[:look_pos]).start()
                seqB = seqB[:l_base] + '-' + seqB[(l_base+1):ins_pos] + base + seqB[(ins_pos+1):]

    return (seqA, seqB, score, begin, end)
## ======================================================================
def shift_to_left_chop(align, matchLen=1):
    ''' Take the result from Bio.pairwise2.align.global**, shift the indels in the homopolymer to the leftmost positions.
    
        Input:  The result (align) is a list of tuples: (seqA, seqB, score, begin, end). seqA and seqB are strings showing the alignment between the sequences. score is the score of the alignment. begin and end are indexes into seqA and seqB that indicate the where the alignment occurs.
                Consider one tuple in the list as input to this function.
                Note: We'll consider that seqA is part of the reference sequence, and seqB is the short read sequence, and the global alignment is done, with no gap penalty at the ends for seqB, but gap penalty at the ends of seqA

        Output: The alignment after shifting
    '''
    #print "align = ", align
    seqA, seqB, score, begin, end = align
    # First find all the insertion positions, ignoring the opening and ending gaps in seqB
    first_non_gap = 0 if seqB[0]!='-' else re.search(r'^[-]+',seqB).end()
    last_non_gap = len(seqB) if seqB[-1]!='-' else re.search(r'[-]+$',seqB).start()

    ins_char = re.compile('-')

    mapping_type_array = zeros(len(seqA),dtype=int32)
    for i in xrange(len(seqA)):
        if seqA[i] == '-' and seqB[i] != '-':
            mapping_type_array[i] = 1 # insertion into seqA
        elif seqB[i] == '-' and seqA[i] != '-':
            mapping_type_array[i] = 2 # deletion from seqA
        elif seqA[i] != seqB[i]:
            mapping_type_array[i] = 3 # mismatch

    match_pos = where(mapping_type_array == 0)[0] # matching positions
    matching_start, matching_end = get_gaps(match_pos) # matching region start and end positions
    matching_lens = matching_end - matching_start # matching regions' lengths
    long_match_ind = where(matching_lens >= matchLen)[0] # long match regions (longer than specified matchLen)
    long_start = matching_start[long_match_ind] # start and end coordinates of the long regions
    long_end = matching_end[long_match_ind]
    shift_align = None
    if long_start[0] != 0:
        seq1 = seqA[:long_start[0]]
        seq2 = seqB[:long_start[0]]
        #print seq1
        #print seq2
        if seq1.count('-') != len(seq1) and seq2.count('-') != len(seq2): # one sequence is all insertion
            seq1 = re.sub('-','',seq1)
            seq2 = re.sub('-','',seq2)
            myalign = pairwise2.align.globalms(seq1, seq2, 0, -1, -0.9, -0.9, penalize_end_gaps=[True, True])
            shift_align = [item for item in pick_align(myalign, False)[0]][:3]
        else:
            shift_align = [seq1, seq2, -0.9*len(seq1)]
        #print format_alignment(*(shift_align + [0, len(seq1)]))

    for i in xrange(len(long_match_ind)-1):
        seq1 = seqA[long_start[i] : long_start[i+1]]
        seq2 = seqB[long_start[i] : long_start[i+1]]
        #print seq1
        #print seq2
        if seq1.count('-') != len(seq1) and seq2.count('-') != len(seq2): # one sequence is all insertion
            seq1 = re.sub('-','',seq1)
            seq2 = re.sub('-','',seq2)
            myalign = pairwise2.align.globalms(seq1, seq2, 0, -1, -0.9, -0.9, penalize_end_gaps=[True, True])
            myalign = [item for item in pick_align(myalign, False)[0]]
        else:
            myalign = [seq1, seq2, -0.9*len(seq1), 0, len(seq1)]
        if shift_align is None:
            shift_align = myalign[:3]
        else:
            shift_align = [shift_align[j] + myalign[j] for j in xrange(3)]
        #print format_alignment(*(shift_align + [0, len(seq1)]))

    if long_end[-1] != len(seqB) - 1: # last base is not a match
        seq1 = seqA[long_start[-1]:len(seqA)]
        seq2 = seqB[long_start[-1]:len(seqB)]
        #print seq1
        #print seq2
        if seq1.count('-') != len(seq1) and seq2.count('-') != len(seq2): # one sequence is all insertion
            seq1 = re.sub('-','',seq1)
            seq2 = re.sub('-','',seq2)
            myalign = pairwise2.align.globalms(seq1, seq2, 0, -1, -0.9, -0.9, penalize_end_gaps=[True, True])
            myalign = [item for item in pick_align(myalign, False)[0]]
        else:
            myalign = [seq1, seq2, -0.9*len(seq1), 0, len(seq1)]
        if shift_align is None:
            shift_align = myalign[:3]
        else:
            shift_align = [shift_align[j] + myalign[j] for j in xrange(3)]
        #print format_alignment(*(shift_align + [0, len(seq1)]))
    else:
        if shift_align is not None:
            shift_align = [shift_align[0] + seqA[long_start[-1]:] , shift_align[1] + seqB[long_start[-1]:], shift_align[2]]
        else:
            shift_align = [seqA[long_start[-1]:] , seqB[long_start[-1]:], 0]

    shift_align += [0]
    shift_align += [len(shift_align[0])]
    return shift_align
## ======================================================================
def pick_align(align_list,trim=True):
    ''' From a list of equivalent alignments between 2 sequences using dynamic progamming (Needleman-Wunch), pick the one whose indel positions are the most left
        Input:  align_list - list of tuples output from Bio.pairwise2.globalXX
                trim - switch to trim the opening and closing gaps in seqB
        Output: align - one of the tuple in the list (seqA, seqB, score, first_non_gap_pos, last_non_gap_pos), last two elements were begin and end originally
    '''
    leftmost_indel_pos = (-1,-1)
    bestalign = ''
    First_non_gap = 0
    Last_non_gap = -1
    #if len(align_list) == 1000:
    #    sys.stderr.write("{} alignments\n".format(len(align_list)))
    for align in align_list:
        seqA, seqB, score, begin, end = align
        # First find all the insertion positions, ignoring the opening and ending gaps in both sequences
        if trim:
            first_non_gap = 0 if seqB[0]!='-' else re.search(r'^[-]+',seqB).end()
            last_non_gap = len(seqB) if seqB[-1]!='-' else re.search(r'[-]+$',seqB).start()

            ## Only look at the aligned part
            seqA = seqA[first_non_gap:last_non_gap]
            seqB = seqB[first_non_gap:last_non_gap]

        ins_char = re.compile('-')
        # insert positions in sequence A and B
        insA = [m.start() for m in ins_char.finditer(seqA)]
        insB = [m.start() for m in ins_char.finditer(seqB)]
        this_indel_pos = (sum(insA+insB), sum(insB)) # use two sums as the measure, one for both sequences and one for the query sequence
        #print "indel positions: ", this_indel_pos
        if leftmost_indel_pos == (-1,-1):
            leftmost_indel_pos = this_indel_pos
            bestalign = (seqA, seqB, score, 0, len(seqA))
            if trim:
                First_non_gap = first_non_gap
                Last_non_gap = last_non_gap
        elif this_indel_pos < leftmost_indel_pos:
            leftmost_indel_pos = this_indel_pos
            bestalign = (seqA, seqB, score, 0, len(seqA))
            if trim:
                First_non_gap = first_non_gap
                Last_non_gap = last_non_gap
    if trim:
        return [bestalign, First_non_gap, Last_non_gap] # return the list of aligns whose indel positions are the leftmost
    else:
        return [bestalign, 0, len(seqA)]
## ======================================================================
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
    if seqB[0] == '-' and seqA[0] == '-': # both deletion, will be ignored in the sequences and CIGAR string
        pass
    elif seqA[0] == seqB[0]: # match
        mode = '='
        cigar_dict['='] = 1
    elif seqA[0] == '-': # insertion
        mode = 'I'
        cigar_dict['I'] = 1
    elif seqB[0] == '-': # deletion
        mode = 'D'
        cigar_dict['D'] = 1
    elif seqA[0] != seqB[0]: # mismatch
        mode = 'X'
        cigar_dict['X'] = 1
    else:
        sys.stderr.write('{} <-> {} in alignment! \n'.format(seqA[0], seqB[0]))

    for i in xrange(1,len(seqA)):
        if seqA[i] == '-' and seqB[i] == '-':
            pass
        elif seqA[i] == seqB[i]:
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
    return cigar_string,first_non_gap,last_non_gap
## ======================================================================
def get_good_regions(ref_bps, rSeq, minGoodLen=1000, minCV=1):
    ''' From mapping results from Illumina reads to a PacBio read, find the regions of PacBio read that satisfy these conditions:
        1.  Every base pair in the region is covered by at least minCV reads
        2.  The contiguous region has length >= minGoodLen
        
        Input:  ref_bps - dictionary for the match/mismatch and the deletion positions
                rSeq - PacBio sequence
                minGoodLen - length threshold for the region contiguously covered by minCV 
                minCV - minimum coverage depth for a base pair to be considered "covered"
        Output: good_regions - list of tuples, each tuple contain the begin and end coordinates of a good region
                cov_bps - number of basepairs covered by the reads
    '''
    alphabet = 'ACGTD' # all the possible base pairs at a position
    rLen = len(rSeq)
    # if the whole sequence is shorter than 1000 bps, do not consider this sequence
    if rLen < minGoodLen:
        return [], [], None

    cov_depths = []
    good_regions = []
    # loop through all the positions of the reference sequence
    for i in xrange(rLen):
        non_ins_bps = ref_bps[i] # non-insertion base calls
        cov_depth = sum(non_ins_bps) # coverage depth at this position
        cov_depths.append(cov_depth) # append the new coverage depth to the vector
    cov_bps = sum([1 for cv in cov_depths if cv > 0]) # number of bases that are covered 
    avg_cov_depth = sum(cov_depths) / float(cov_bps) if cov_bps != 0 else 0 # average coverage depth for the covered bases
    low_CV_pos = [-1] + [ i for i in xrange(rLen) if cov_depths[i] < minCV ] + [rLen] # indices of the lower coverage bases (coverage depth lower than the specified minCV)
    for i in xrange(1, len(low_CV_pos)):
        if low_CV_pos[i] - low_CV_pos[i-1] >= minGoodLen:
            good_regions.append((low_CV_pos[i-1]+1, low_CV_pos[i])) # found a good region (begin, end) where rSeq[begin:end] is long enough and covered well
    return good_regions, cov_bps, avg_cov_depth
## ======================================================================
def get_poly_pos(ref_bps, ref_ins_dict, region=None):
    ''' Get the consensus and polymorphic positions where the insertion happened, in the specified region, return as dictionaries
        Input:  ref_bps, ref_ins_dict - output objects from read_and_process_sam, noninsertion and insertion information for all the positions from the Illumina reads
                region - (begin, end) tuple of a region to consider
        Output: poly_bps - polymorphic non-insertion positions
                poly_ins - polymorphic insertion positions
                consensus_bps - nonpolymorphic non-insertion positions, get the consensus
                consensus_ins - nonpolymorphic insertion positions, get the consensus
                * note: all the above four objects are dictionaries. Dictionaries for the insertion positions use tuple (insert position, index of inserted base) as key, since list cannot be used as dictionary key.
                cvs - coverage depths across all the base pairs
    '''
    if region is None: # region to consider, 0-based index
        region = (0, len(ref_bps))
    poly_bps = dict() # list of tuples
    poly_ins = dict()
    consensus_bps = dict()
    consensus_ins = dict()
    alphabet = 'ACGTD' # all the possible base pairs at a position
    cvs = [] # list of coverage depths

    # consider all the positions in the specified region
    for pos in xrange(region[0], region[1]):
        # check the match/mismatch/deletion first
        cv = sum(ref_bps[pos]) # coverage depth
        cvs.append(cv) # vector of coverage depths, for each position in the region
        base_calls = [ alphabet[i] for i in xrange(5) if ref_bps[pos][i] >= minCover(cv) ]
        if len(base_calls) > 1: # more than 1 base call with enough read support
            poly_bps[pos] =  base_calls
        elif len(base_calls) == 1:
            consensus_bps[pos] = base_calls[0]
        else: # no base has more than required number of coverage, treat base call with 1 read support as valid
            base_calls = [ alphabet[i] for i in xrange(5) if ref_bps[pos][i] > 0 ]
            if len(base_calls) == 1: # if they all agree (low CV region) then use the consensus base
                consensus_bps[pos]= base_calls[0]
            else: # if they do not agree, treat each one as possible/valid. because all of them should have very low read support now passing 
                poly_bps[pos] = base_calls
                #consensus_bps.append((pos, alphabet[ref_bps[pos].index(max(ref_bps[pos]))]))

        # check the insertion positions now
        if pos in ref_ins_dict: 
            for i in xrange(len(ref_ins_dict[pos])): # look at all the bases inserted at this position
                pos_list = ref_ins_dict[pos][i]
                base_calls = [ alphabet[j] for j in xrange(4) if pos_list[j] >= minCover(cv) ]
                if len(base_calls) > 1: # polymorphic insertion position
                    poly_ins[ (pos,i) ] =  base_calls
                elif len(base_calls) == 1: # consensus insertion position
                    consensus_ins[ (pos,i) ] =  base_calls[0]
                # if none of them passes the threshold, the insertion position won't be considered existent
    return poly_bps, poly_ins, consensus_bps, consensus_ins, cvs
## ======================================================================
def ref_extension(poly_bps, poly_ins, consensus_bps, consensus_ins, rseq, region=None, print_width = 100, verbose=False):
    ''' Extend the reference sequence, insert the insertion positions in the reference sequence, and return correspondence between old positions and their new positions in the extended sequence
        
        Input:  poly_bps, poly_ins, consensus_bps, consensus_ins - output from previous function (get_poly_pos)
                rseq - reference sequence
                region - (begin,end) of the region to work on

        Output: newSeq - new extended reference sequence corresponding to the region
                bp_pos_dict - dictionary old position => new position in the extended sequence
                ins_pos_dict - same dictionary for the insertion positions
    '''
    alphabet = 'ACGTD' # all the possible base pairs at a position
    rLen = len(rseq)
    if region is None: # if region is not specified, get it from poly_bps and consensus_bps, from the leftmost covered position to the rightmost covered position
        begin = min(consensus_bps.keys() + poly_bps.keys())
        end = max(consensus_bps.keys() + poly_bps.keys())
        region = (begin,end+1)
    else:
        begin = region[0]
        end = region[1] - 1

    ins_pos = sorted(consensus_ins.keys() + poly_ins.keys()) # all the insertion positions, sorted list of tuples
    bp_pos = sorted(consensus_bps.keys() + poly_bps.keys()) # all the non-insertion positions, sorted list of integers

    # initilalization
    ins_pos_dict = dict()
    bp_pos_dict = dict()
    inserted_bases = 0
    current_pos = begin

    newSeq = '' # new sequence with consensus base call non-insertion positions updated, and dashes for insertions and deletions
    oldSeq = '' # old sequence with just dashes for insertions
    matching_string = '' # matching string to print between two sequences for better visualization

    for ins in ins_pos: # loop through all the insertion positions, these are the places changing the length of the sequence
        for mypos in xrange(current_pos, ins[0]): # non-insertion positions
            if mypos not in consensus_bps: # position has polymorphism or doesn't have enough short read coverage, keep the original base call
                newSeq += rseq[mypos]
                oldSeq += rseq[mypos]
                matching_string += '^'
            else: # consensus position
                newSeq += consensus_bps[mypos] # new sequence takes the consensus position
                oldSeq += rseq[mypos]
                if newSeq[-1] == oldSeq[-1]:
                    matching_string += '|'
                else:
                    matching_string += 'X'
            bp_pos_dict[mypos] = mypos + inserted_bases # old position => new position in the extended ref sequence
            #print "bps_dict: {} => {}".format(mypos, mypos+inserted_bases)
        # now look at the insertion position
        if ins in consensus_ins:
            newSeq += consensus_ins[ins] # insert bases
            matching_string += ' '
        else:
            newSeq += poly_ins[ins][0] # insert bases
            matching_string += '^'
        oldSeq += '-'
        ins_pos_dict[ins] = ins[0] + inserted_bases
        #print "ins_dict: ({},{}) => {}".format(ins[0], ins[1], ins[0]+inserted_bases)
        inserted_bases += 1 # update total number of inserted bases so far
        current_pos = ins[0]

    # positions after the last insertion
    for mypos in xrange(current_pos, region[1]):
        if mypos not in consensus_bps: # position has polymorphism or doesn't have enough short read coverage
            newSeq += rseq[mypos]
            oldSeq += rseq[mypos]
            matching_string += '^'
        else:
            newSeq += consensus_bps[mypos] # if there is consensus
            oldSeq += rseq[mypos]
            if newSeq[-1] == oldSeq[-1]:
                matching_string += '|'
            else:
                matching_string += 'X'
        bp_pos_dict[mypos] = mypos + inserted_bases # old position => new position in the extended ref sequence
        #print "bps_dict: {} => {}".format(mypos, mypos+inserted_bases)

    # print the alignment between the original sequence and the newly proposed sequence
    # for the two different sequences, print the coordinates on each row, for both of them.
    if verbose:
        sys.stdout.write("Extend the PacBio sequence with the insertion positions in the specified region with extended sequence below original sequence.\n")
        sys.stdout.write("X means consensus base call different from original\n")
        sys.stdout.write("^ means positions is polymorphic\n\n")
        for i in xrange(len(newSeq) / print_width + 1):
            sys.stdout.write('{:>5}\t{}\t{:<5}\n'.format(begin + i*print_width - oldSeq[:i*print_width].count('-'), oldSeq[i*print_width:(i+1)*print_width], min(begin+(i+1)*print_width - oldSeq[:(i+1)*print_width].count('-'),len(oldSeq) + begin -oldSeq.count('-'))))
            sys.stdout.write('{:>5}\t{}\t{:<5}\n'.format(' ', matching_string[i*print_width:(i+1)*print_width],' ' ))
            sys.stdout.write('{:>5}\t{}\t{:<5}\n\n'.format(begin + i*print_width, re.sub('D','-',newSeq[i*print_width:(i+1)*print_width]), min(begin+(i+1)*print_width,len(newSeq) + begin) ))

    return newSeq, bp_pos_dict, ins_pos_dict
## ======================================================================
def update_pos_info(poly_bps, poly_ins, consensus_bps, consensus_ins, bp_pos_dict, ins_pos_dict):
    ''' Given the correspondence between the positions in the original sequence and the positions in the extended sequence, update the polymorphic and consensus information for both insertion and non-insertion positions.
        Input:  poly_bps, poly_ins, consensus_bps, consensus_ins, bp_pos_dict, ins_pos_dict
        Output: poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext
    '''
    poly_bps_ext = dict()
    poly_ins_ext = dict()
    consensus_bps_ext = dict()
    consensus_ins_ext = dict()
    for i in poly_bps:
        poly_bps_ext[ bp_pos_dict[i] ] = poly_bps[i]
    for i in poly_ins:
        poly_ins_ext[ ins_pos_dict[i] ] = poly_ins[i]
    for i in consensus_bps:
        consensus_bps_ext[ bp_pos_dict[i] ] = consensus_bps[i]
    for i in consensus_ins:
        consensus_ins_ext[ ins_pos_dict[i] ] = consensus_ins[i]
    return poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext
## ======================================================================
def make_type_array(poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext,verbose=False):
    ''' Construct an array that indicates the type each position. 0: consensus non-insertion; 1: consensus insertion; 2: polymorphic non-insertion; 3: polymorphic insertion.
        Input:  position information list (including poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext) for the extended sequence
        Output: an array with 0,1,2,3 as entries and with length equal to the good region
                (start, end) - start and end of the region in extended sequence
    '''
    poly_bps_pos = array( [i for i in poly_bps_ext], dtype=int32 ) # extract positions in the big sequence, for each type
    consensus_bps_pos = array( [i for i in consensus_bps_ext], dtype=int32 )
    poly_ins_pos = array( [ i for i in poly_ins_ext], dtype=int32 )
    consensus_ins_pos = array( [ i for i in consensus_ins_ext], dtype=int32 )
    all_pos = hstack((poly_bps_pos, consensus_bps_pos, poly_ins_pos, consensus_ins_pos)) # all the positions

    start = min(all_pos) # start and end of the region currently covered by the input
    end = max(all_pos) + 1
    if verbose:
        sys.stdout.write('Making type array for the new sequence from {} to {}.\n'.format(start, end))
    type_array = zeros(end - start, dtype=int32) # initialize all to type 0 
    if len(consensus_ins_pos) > 0: # fill in different types
        type_array[consensus_ins_pos - start] = 1
    if len(poly_bps_pos) > 0:
        type_array[poly_bps_pos - start] = 2
    if len(poly_ins_pos) > 0:
        type_array[poly_ins_pos - start] = 3
    return type_array, (start, end)
## ======================================================================
def make_ref_array(consensus_bps_ext, consensus_ins_ext, type_array,ext_region):
    ''' Construct array for the reference sequence, with consensus positions' entries filled in.
        Polymorphic positions are filled with 0s now.
        Input:  consensus information (consensus_bps_ext, consensus_ins_ext)
                type_array - array including type of each position
                ext_region - (start, end) of the region on the extended sequence
        Output: ref_array - array with size 5 times size of type_array, each position has 5 options (ACGTD)
    '''
    alphabet = 'ACGTD'
    ref_array = zeros( len(type_array) * 5, dtype = int32 )
    start = ext_region[0]
    for i in where(type_array == 0)[0]:
        base = [ j[1] for j in consensus_bps_ext if j[0]-start == i][0]
        ref_array[ i*5 + alphabet.index(base) ] = 1
    for i in where(type_array == 1)[0]:
        base = [ j[1] for j in consensus_ins_ext if j[0]-start == i][0]
        ref_array[ i*5 + alphabet.index(base) ] = 1
    return ref_array
## ======================================================================
def array_to_seq(seq_array, miss_char='.'):
    ''' Convert 0-1 1d array back to nucleotide sequence.
        Input:  seq_array - 1d array with 0 and 1 as entries. length has to be multiple of 5, there can only be one 1 every 5 positions from start
                miss_char - character for the bases not covered by this read (all 5 0s for this position)
        Output: seq_short - corresponding DNA sequence with deletion removed
                seq_long - corresponding DNA sequence with deletion replaced by '-'
    '''
    if len(seq_array) % 5 != 0: # array length has to be a multiple of 5
        sys.exit("length of the array ({}) is not a multiple of 5!! \n".format(len(seq_array)))
    else:
        seqLen = len(seq_array) / 5 # sequence length
        seq_array = seq_array.reshape(-1,5) # reshape to an array with 5 columns
        seq = [miss_char] * seqLen # initialize the sequence to return, as a list
        for i in xrange(5): # for each base (ACGTD) look for their positions
            positions = where(seq_array[:,i] == 1)[0]
            seq = [ alphabet[i] if x in positions else seq[x] for x in xrange(seqLen) ] # update the sequence with the specified base call
        seq = ''.join(seq) # join list to a string
        seq_short = re.sub('D','',seq) # remove the deletion positions from sequence
        seq_long = re.sub('D','-',seq)
        return seq_short, seq_long
## ======================================================================
def print_seqs(seq1, seq2, print_width = 100):
    ''' Print two sequences next to each other, indicate the part that has changed with lower case DNA letters.
        Input:  seq1, seq2 - two sequences to print/compare/display
                print_width - number of nucleotides to print each line
        Output: standard output of the two sequences next to each other
    '''
    if len(seq1) != len(seq2):
        sys.exit("Input sequences do not have same length! \n")
    matching_string = ''
    for i in xrange(len(seq1)):
        if seq1[i].upper() == seq2[i].upper():
            matching_string += '|'
        else:
            matching_string += 'X'
            #seq1 = seq1[:i] + seq1[i].lower() + seq1[i+1:]
            #seq2 = seq2[:i] + seq2[i].lower() + seq2[i+1:]

    for i in xrange(len(seq1) / print_width + 1):
        sys.stdout.write('{:>5}\t{}\t{:<5}\n'.format(i*print_width, seq1[i*print_width:(i+1)*print_width], min((i+1)*print_width,len(seq1) )))
        sys.stdout.write('{:>5}\t{}\t{:<5}\n'.format(' ', matching_string[i*print_width:(i+1)*print_width],' ' ))
        sys.stdout.write('{:>5}\t{}\t{:<5}\n\n'.format(i*print_width - seq2[:(i*print_width)].count('-'), seq2[i*print_width:(i+1)*print_width], min((i+1)*print_width,len(seq1) ) - seq2[:min((i+1)*print_width,len(seq1) )].count('-') ))
## ======================================================================
def make_read_array1d(read_string, bp_pos_dict, ins_pos_dict, type_array, poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext, ext_region=None):
    ''' Make 1d array for a particular read from its string (key of dictionary readinfo).
        Incorporate the base calling from other reads at the same position.
        Considering how many bases a read need to change to comply with the consensus calls, if too many, keep its own base calls (or discard it)

        Input:  read_string - read string (key of dictionary readinfo)
                bp_pos_dict, ins_pos_dict - correspondence between positions from original to extended ref sequence
                length - length of the reference sequence ( good region )
                poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext - base calling information for all the positions
                ext_region - (start, end) of the region on the extended sequence

        Output: read_array1d - 1d array for this particular read string
    '''
    if ext_region is None:
        all_pos = bp_pos_dict.values() + ins_pos_dict.values() # values instead of keys, since it's the "extended" region
        ext_region = (min(all_pos), max(all_pos)+1 ) # the end is one bigger than the 0-based end
    start, end = ext_region

    read_array1d = zeros( len(type_array) * 5, dtype=int32 ) # initialization
    bpstring = read_string.split(":")[0] # non-insertion base calls
    insstring = read_string.split(":")[1] # insertion position's base calls

    # all the positions of a read that lay in the specified region (covered by the dicts)
    bp_positions = []
    ins_positions = []
    for pos in map(int, re.findall('\d+',bpstring)):
        if pos in bp_pos_dict:
            bp_positions.append(bp_pos_dict[pos])
    ins_chars = re.findall('\D+',insstring)
    ins_pos = re.findall('\d+',insstring)
    for i in xrange(len(ins_pos)):
        for j in xrange(len(ins_chars)):
            if (ins_pos[i], j) in ins_pos_dict:
                ins_positions.append(ins_pos_dict[(ins_pos[i], j)])
    map_positions = bp_positions + ins_positions
    #print len(map_positions), " covered positions: ", map_positions

    if len(map_positions) == 0: # the read does not cover any (non-insertion and insertion) of the positions in the region
        return read_array1d
    else:
        start_pos = min(map_positions) # start and end positions of the covered region for this read
        end_pos = max(map_positions)

        changed_bases = 0

        # If this read covers some non-insertion positions, do the following:
        if len(bpstring) > 0:
            positions = map(int, re.findall('\d+',bpstring)) # positions
            bases = re.findall('\D+',bpstring) # base calling
            for i in xrange(len(positions)):
                if positions[i] in bp_pos_dict:
                    position = bp_pos_dict[positions[i]] # position in the extended sequence
                    ## note: type_array's index starts with 0, instead of "start", so does read_array1d
                    relative_pos = position - start
                    if type_array[relative_pos] == 0: # consensus non-inseriton position, read's call should be the same as the consensus call
                        c_bp = consensus_bps_ext[position] # consensus base
                        read_array1d[ (relative_pos)*5 + alphabet.index(c_bp) ] = 1
                        if bases[i] != c_bp:
                            #sys.stdout.write("{} is not {}, changed_bases plus 1\n".format(bases[i], c_bp))
                            changed_bases += 1
                    if type_array[relative_pos] == 2: # polymorphic non-insertion position, check if read's call is one of the possible calls
                        p_bps = poly_bps_ext[position] # polymorphic base
                        if bases[i] in p_bps: # if read's call is one of the possible calls, do not change it
                            read_array1d[ (position-start)*5 + alphabet.index(bases[i]) ] = 1
                        else: # if read's call is not among the possible calls, it's considered as an error, so it could be either of the possible calls
                            for p_bp in p_bps:
                                read_array1d[ (position-start)*5 + alphabet.index(p_bp) ] = 1
                            #sys.stdout.write("{} is not in {}, changed_bases plus 1\n".format(bases[i], ''.join(p_bps)))
                            changed_bases += 1
                            
        # If this read covers some insertion positions, do the following:            
        if len(insstring) > 0:
            positions = map(int, re.findall('\d+',insstring))
            ins_positions = [] # insert positions, starting from 0
            bases = re.findall('\D+',insstring)
            for i in xrange(len(positions)): # for each insertion position in this read
                for j in xrange(len(bases[i])):
                    if (positions[i],j) in ins_pos_dict:
                        ins_position = ins_pos_dict[(positions[i],j)]
                        ins_positions.append(ins_position - start) # all the insertion positions the read covers
                        if (ins_position - start) < len(type_array): # in some cases, when there is more than 1 base pair inserted, the index could go out of range
                            if type_array[ins_position - start] == 1: # consensus insertion position
                                c_bp = consensus_ins_ext[ins_position]
                                read_array1d[ (ins_position-start)*5 + alphabet.index(c_bp) ] = 1
                                if bases[i][j] != c_bp:
                                    #sys.stdout.write("{} is not {}, changed_bases plus 1\n".format(bases[i][j], c_bp))
                                    changed_bases += 1

                            if type_array[ins_position - start] == 3: # polymorphic insertion position
                                p_bps = poly_ins_ext[ins_position]
                                if bases[i][j] in p_bps:
                                    read_array1d[ (ins_position-start)*5 + alphabet.index(bases[i][j]) ] = 1
                                else:
                                    for p_bp in p_bps:
                                        read_array1d[ (ins_position-start)*5 + alphabet.index(p_bp) ] = 1
                                    changed_bases += 1
        #print changed_bases

        if changed_bases >= 5: # if this read has at least 5 disagreements with the consensus/polymorphic calls, discard it (set the coverage to nothing)
            return zeros( len(type_array) * 5, dtype=int32 )
        else:
            # check if this read missed any insertion position
            for ins_pos in hstack((where(type_array==3)[0], where(type_array==1)[0])):
                if ins_pos >= start_pos-start and ins_pos <= end_pos-start and ins_pos not in ins_positions: # if there is no insertion for this read at the insertion position, put a deletion
                    read_array1d[ins_pos * 5 + 4 ] = 1
            return read_array1d
## ======================================================================
def make_read_array(readinfo, bp_pos_dict, ins_pos_dict, type_array,  poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext, ext_region=None):
    ''' Make 2d array for all reads from readinfo dictionary, dim1 of the array is equal to the length of the dictionary
        Input:  readinfo - dictionary (read string => count of the same read string)
                bp_pos_dict, ... - same as function make_read_array1d's input for each read string
        Output: read_array - 2d array that include all reads' base call information
                read_counts - 1d array that stores the counts of reads corresponding to each row of read_array
    '''
    read_array = zeros( (len(readinfo), len(type_array)*5), dtype = int32 ) # initialize the 2d array to return
    read_counts = zeros( len(readinfo), dtype = int32)
    i = 0
    for read_string in sorted(readinfo): # the rows in the array are ordered by the corresponding read_string
        read_array[i,:] = make_read_array1d(read_string, bp_pos_dict, ins_pos_dict, type_array, poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext, ext_region)
        read_counts[i] = len(readinfo[read_string])
        i += 1
    read_array = read_array[where(read_array.sum(axis=1) !=0)[0],:] # only the reads that cover some bases of the region
    base_count = read_array.sum(axis=0).reshape(-1,5)
    covered_bases = where(base_count.sum(axis=1) != 0)[0]
    return read_array[where(read_array.sum(axis=1) !=0)[0],amin(covered_bases)*5:amax(covered_bases+1)*5], read_counts
## ======================================================================
def is_compatible(array1, array2):
    ''' Find out if the base calls of 2 reads at a certain position are compatible or not, given the 0-1 vectors representing the base calls. For now, the length of the vectors is the same (5 for ACGTD).
        Input:  array1, array2 - length 5 vectors, each corresponds to the base call of a read at one position
        Output: boolean value - true if they are compatible, false if not
    '''
    # only 5 entries each array in array1 and array2
    same_call_pos = where(bitwise_and( array1, array2) == 1)[0] # find number of nucleotides that both reads call
    if len(same_call_pos) > 0: # if they agree on at least 1 nucleotide, they are compatible, otherwise they are not
        return True
    else:
        return False
## ======================================================================
def are_reads_compatible(read_array1d1, read_array1d2):
    ''' Find if two reads are compatible, given their 1d array with 0 and 1 entries.
        Input:  read_array1d1, read_array1d2 - 0-1 vectors for 2 reads
        Output: boolean value - true if they are compatible, false if not
    '''
    r1 = read_array1d1.reshape(-1,5) # change the vector to 2-d arrays
    r2 = read_array1d2.reshape(-1,5)
    overlap_pos = intersect1d( where(r1==1)[0], where(r2==1)[0] ) # find the overlap positions of the 2 reads
    # if 2 reads do not overlap then they are compatible
    if len(overlap_pos) == 0:
        return True
    else: # otherwise, check position by position and see if any is not compatible
        for pos in overlap_pos:
            if not is_compatible(r1[pos,:], r2[pos,:]): # if non-compatible position is found, return False immediately and quit function
                return False
    return True
## ======================================================================
def compatible_mat(read_array):
    ''' Construct (upper triangular) compatibility matrix for pairwise compatibility of the reads whose info is stored in the read_array.
        Input:  read_array- 2d array from all the reads 
        Output: Cmat - compatibility array for all the reads in the 2d array
    '''
    nread = read_array.shape[0] # number of rows, which is number of the reads 
    Cmat = zeros((nread, nread),dtype=int32) # initialize to be all 0 - incompatible
    for i in xrange(nread):
        for j in xrange(i+1, nread):
            if are_reads_compatible(read_array[i],read_array[j]): # set the entry for compatible pairs to 1
                Cmat[i,j] = 1
    return array(Cmat) # conver to array
## ======================================================================
def cov_bps(read_array1d):
    ''' Find number of bases covered by an array 
        Input: read_array1d - 1 dimensional array for a read
        Output: len(cov_pos) - number of covered base pairs by this read
    '''
    r = read_array1d.reshape(-1,5) # convert to array with 5 columns
    cov_pos = unique( where( r==1 )[0] ) # positions covered by this read, remove repeat when a read gives ambiguous base calls
    return len(cov_pos) # number of base pairs covered by the given read (array)
## ======================================================================
def is_read_compatible(ref_array, read_array1d):
    ''' Given a reference array and an array for a read, determine if they are compatible or not.
        Input:  ref_array - vector of 0-1 for base calls on the reference sequence
                read_array1d - 1d array from one read
        Output: True if read is compatible, False otherwise
    '''
    r = read_array1d.reshape(-1,5)
    cov_pos = unique( where( r.sum(axis=1)>0 )[0] ) # positions covered by this read
    #print "number of covered bases is ", len(cov_pos) # DEBUG
    if dot(ref_array, read_array1d) == len(cov_pos): # there should only be one 1 for the length 5 vector of a base pair position
        return True
    else:
        return False
## ======================================================================
def get_compatible_reads(ref_array, read_array):
    ''' Given a candidate reference sequence, find reads that are compatible with this candidate
        Input:  ref_array - vector of 0-1 for base calls on the reference sequence
                read_array - 2d array from all the reads
        Output: compatible_ind - array of indices for compatible reads, corresponding to the read_array
    '''
    compatible_ind = []
    for i in xrange(read_array.shape[0]): # check every read
        if sum(read_array[i,:]) > 0 and is_read_compatible(ref_array, read_array[i,:]): # ignore reads that have all 0 vector
            compatible_ind.append(i)
    return array(compatible_ind, dtype=int32)
## ======================================================================
def get_overlapLen(ref_array, read_array, Cvec=None):
    ''' Given read array, current reference array, and the indices of compatible reads, get the overlap length matrix.
         Input: ref_array - array for the currently proposed PacBio sequence
                read_array - array with read information
                Cvec - indices for the compatible reads with the PacBio sequence, doesn't have to be given
        Output: overlap_mat - matrix with overlap lengths
    '''
    if Cvec is None:
        Cvec = get_compatible_reads(ref_array, read_array)
    if len(Cvec) == 0:
        return None
    nreads = len(Cvec)
    new_array = read_array[ Cvec, : ] # get only the compatible reads
    overlap_mat = zeros(shape=(nreads,nreads),dtype=int32) # initialize the overlap length matrix, with each compatible read having one row and one column
    zero_columns = where(ref_array==0)[0] # force only 1 base call for the reads that have ambiguous (or error) base calls, and they agree with the reference's base call
    new_array[ : , zero_columns] = 0
    start_pos_vec = zeros(nreads,dtype=int32)
    end_pos_vec = zeros(nreads,dtype=int32)
    for i in xrange(nreads): # starting position of the read
        start_pos_vec[i] = where(new_array[i,:] > 0)[0][0]/5
        end_pos_vec[i] = where(new_array[i,:] > 0)[0][-1]/5
    for i in xrange(nreads - 1):
        for j in xrange(i + 1, nreads):
            overlap_len = dot(new_array[i,:], new_array[j,:])
            if start_pos_vec[i] < start_pos_vec[j] and end_pos_vec[i] < end_pos_vec[j]: # positive if the first read's starting position is smaller than the second read's starting position
                overlap_mat[i,j] = overlap_len
                overlap_mat[j,i] = - overlap_mat[i,j]
            elif start_pos_vec[i] > start_pos_vec[j] and end_pos_vec[i] > end_pos_vec[j]:
                overlap_mat[i,j] = -overlap_len
                overlap_mat[j,i] = - overlap_mat[i,j]
    return overlap_mat
## ======================================================================
def cov_vec(read_array1d):
    return apply_along_axis(max, 1, read_array1d.reshape(-1,5))
## ======================================================================
def get_overlapMat(read_array):
    ''' Given read array, get the overlap length matrix between all the reads. If two reads are not compatible, their overlap length will be set to 0
        Too slow...
        Input: read_array - array with read information
        Output: overlap_Mat - matrix with overlap lengths between all reads
    '''
    nreads = read_array.shape[0] # number of reads in the read_array
    overlap_Mat = zeros(shape=(nreads,nreads),dtype=int32) # initialize the overlap length matrix, with each read having one row and one column
    cov_mat = apply_along_axis(cov_vec, 1, read_array) # coverage matrix

    # reads that cover at least 1 bp
    cov_reads = where(sum(read_array,1) != 0)[0]
    # start and end positions of each read
    start_pos_vec = zeros(len(cov_reads),dtype=int32)
    end_pos_vec = zeros(len(cov_reads),dtype=int32)
    for i in xrange(len(cov_reads)): # starting and ending positions of the read
        start_pos_vec[i] = where(cov_mat[cov_reads[i],:] > 0)[0][0]
        end_pos_vec[i] = where(cov_mat[cov_reads[i],:] > 0)[0][-1]
    start_sort = argsort(start_pos_vec) # indices when sorted by start of covered positions
    # put the corresponding vectors in this order
    start_pos_vec = start_pos_vec[start_sort]
    end_pos_vec = end_pos_vec[start_sort]
    cov_reads = cov_reads[start_sort]
    for i in xrange(len(start_sort) - 1):
        for j in xrange( i+1, len(start_sort)):
            overlap_len = sum(logical_and(cov_mat[cov_reads[i],], cov_mat[cov_reads[j],])) # overlap length
            agree_len = len(unique(where(logical_and(cov_mat[cov_reads[i],], cov_mat[cov_reads[j],]))[0])) # number of bps where they agree
            if overlap_len == 0:
                continue
            else:
                if end_pos_vec[j] > start_pos_vec[i] and overlap_len == agree_len: # only when one is not included in the other and they are compatible
                    overlap_Mat[cov_reads[i],cov_reads[j]] = overlap_len
                    overlap_Mat[cov_reads[j],cov_reads[i]] = -overlap_len
    #sys.stdout.write("average overlap length is {}\n".format(mean(abs(overlap_Mat)[where(abs(overlap_Mat)>0)]))) ## DEBUG, for soft minOverlap cutoff shortly
    return overlap_Mat
## ======================================================================
def find_maxOverlap(overlap_mat):
    ''' From the overlap length matrix, find the maximum overlap for each read at both left and right directions.
        Input:  overlap_mat - overlap length matrix from get_overlapLen function
        Output: maxOverlap_mat - maximum overlap length array, one row for each read in the overlap_mat
    '''
    maxOverlap_mat = zeros(shape=(overlap_mat.shape[0], 2), dtype = int32) # maximum overlap length for each read, one with reads to its left, one with reads to its right
    ## for each read, do following
    for i in xrange(overlap_mat.shape[0]):
        i_array = hstack((overlap_mat[i,:i],overlap_mat[i, (i+1):])) # remove the diagonal element (overlap with itself)
        #print i_array
        ## LEFT
        left_overlap_array = i_array[ where(i_array < 0)[0]] # negative overlap lengths (with reads to the left of the current read)
        if len(left_overlap_array) > 0: # there are reads to the left that overlap with current read
            left_maxOverlap = amax(abs(left_overlap_array))
        else: # otherwise
            left_maxOverlap = 0

        ## RIGHT
        right_overlap_array = i_array[ where(i_array > 0)[0]] # negative overlap lengths (with reads to the right of the current read)
        if len(right_overlap_array) > 0:# there are reads to the right that overlap with current read
            right_maxOverlap = amax(right_overlap_array)
        else: # otherwise
            right_maxOverlap = 0
        #print left_maxOverlap, "\t", right_maxOverlap
        maxOverlap_mat[i,] = array([left_maxOverlap, right_maxOverlap])

    return maxOverlap_mat
## ======================================================================
def get_reads_name(readinfo, compatible_ind):
    ''' Find the names of reads that are compatible with a given PacBio sequence, given the readinfo dictionary and the compatible array row indices.
        Input:  readinfo - dictionary with read information read_string => list of read names that have this read_string
                compatible_ind - row indices of the array compatible with the PacBio sequence
        Output: reads_name_list - list of names of reads corresponding to the compatible_ind
    '''
    if len(compatible_ind) == 0: # no compatible reads
        return []
    else:
        reads_name_list = [] # initialize an empty list
        info_keys = sorted(readinfo.keys()) # sort the keys of readinfo
        for ind in compatible_ind:
            reads_name_list += readinfo[ info_keys[ind] ]
        return reads_name_list
## ======================================================================
def gap_pos(ref_array, read_array, compatible_ind, minOverlap = 10, minOverlapRatio=0.1):
    ''' Given a reference array and indices (of array) of compatible reads, find the gap positions, i.e. positions that are not covered by any read
        Input:  ref_array - vector of 0-1 for base calls on the reference sequence
                read_array - 2d array from all the reads
                compatible_ind - indices of compatible reads, corresponding to the read_array
        Output: gap_vec - vector with positions that have no coverage from the reads
    '''
    base_pos = where( ref_array == 1)[0] # coordinates corresponding to the bases called by the reference
    ## TODO ## need to exclude the deleted bases from the base positions since they won't exist in the resulted reference sequence
    if len(compatible_ind) == 0: # no compatible reads, all positions are gap positions
        return arange(ref_array.shape[0]/5)
    else:
        sub_read_array = copy(read_array[compatible_ind,:]) # read array with only the compatible reads
        cvec = arange(len(compatible_ind),dtype=int32)
        if minOverlap != -1 or minOverlapRatio != 0: # check minOverlap if at least one of the minOverlap and minOverlapRatio is specified
            iter_count = 1
            trimmed_reads = zeros(shape=(len(compatible_ind),2),dtype=int32)
            while True: # iteratively change the read array if necessary
                overlap_mat = get_overlapLen(ref_array, sub_read_array, cvec)
                if amax(abs(overlap_mat)) == 0:
                    avgOverlap = 0
                else:
                    avgOverlap = mean(abs(overlap_mat)[where(abs(overlap_mat)>0)]) # average overlap length
                cutOverlap = max(minOverlap, avgOverlap * minOverlapRatio) # hard cutoff and soft ratio cutoff, whichever is larger will be used
                #print "iteration ", iter_count
                #print "minimum overlap length: ", cutOverlap

                maxOverlap_mat = find_maxOverlap(overlap_mat)
                small_ind = where( maxOverlap_mat <= cutOverlap) # check if any read need clipping because of small overlap length
                #print small_ind
                #print trimmed_reads
                #print "trim status: ", trimmed_reads[small_ind]
                if len(small_ind[0]) == 0 or all(trimmed_reads[small_ind] == 1):
                #if len(small_ind[0]) == 0 or iter_count > 4:
                    #print "all problematic reads trimmed already"
                    break
                else: # some reads have problem
                    iter_count += 1
                    for i in xrange(len(small_ind[0])):
                        #print i, "==>", small_ind[0][i], ",", small_ind[1][i]
                        if small_ind[1][i] == 0 and trimmed_reads[small_ind[0][i],0] == 0: # left overlap problem
                            trimmed_reads[small_ind[0][i],0] = 1 # mark that this read has been trimmed for the overlap checking
                            if len(where(sub_read_array[small_ind[0][i],:] != 0)[0]) > 0: # if the positions where there are base calls in this read is greater than 0
                                first_pos = where(sub_read_array[small_ind[0][i],:] != 0)[0][0]/5 # first covered position of this read
                                if first_pos > 0 : # if the read is not mapped to the beginning of the region, decrease its coverage accordingly
                                    length = maxOverlap_mat[small_ind[0][i], small_ind[1][i]]
                                    #sys.stdout.write("left: index {} has maximum overlap length {}\n".format(i, length)) #DEBUG
                                    sub_read_array[small_ind[0][i], first_pos*5 : (first_pos + length + 1)*5] = 0 

                        if small_ind[1][i] == 1 and trimmed_reads[small_ind[0][i],1] == 0: # right overlap problem
                            trimmed_reads[small_ind[0][i],1] = 1
                            if len(where(sub_read_array[small_ind[0][i],:] != 0)[0]) > 0:
                                last_pos = where(sub_read_array[small_ind[0][i],:] != 0)[0][-1]/5 # last covered position of this read
                                if last_pos < (len(ref_array)/5 - 1) : # if the read is not mapped to the end of the region, decrease its coverage accordingly
                                    length = maxOverlap_mat[small_ind[0][i], small_ind[1][i]]
                                    #sys.stdout.write("right: index {} has maximum overlap length {}\n".format(i, length))
                                    sub_read_array[small_ind[0][i],(last_pos - length)*5 : (last_pos + 1)*5] = 0 
        base_cov = sub_read_array[ : , base_pos]
        base_cov = base_cov.sum(axis=0) # pick the compatible rows and the correct columns, sum over the columns
        #print "gap positions: ", where(base_cov == 0)[0]
        return where(base_cov == 0)[0]
## ======================================================================
def get_consensus_from_array(read_array):
    ''' get initial starting point for the greedy algorithm, i.e., the consensus sequence found from the read array.
        Input:  read_array - 2d array from all the reads
        Output: ref_consensus - consensus 0-1 array for the reference
    '''
    pos_sum = read_array.sum(axis=0) # take the column sum
    pos_sum_array = pos_sum.reshape(-1,5) # convert to 2d array, with 5 columns
    pos_max = argmax(pos_sum_array, axis=1) # find the index of the base with maximum coverage
    consensus_array = zeros(pos_sum_array.shape,dtype = int32) # initialize the consensus_array to all 0
    consensus_array[ arange(pos_sum_array.shape[0]), pos_max ] = 1 # fill the consensus base coordinates with 1
    return consensus_array.flatten() # flatten 2d array to 1d array
## ======================================================================
def get_gaps(gaps):
    ''' Find the start position of the widest gap. Note: This position itself might not be a polymorphic position to change base.
        Input:  gaps - gap positions (after a ref seq candidate is proposed and compatible reads are extracted, output from function gap_pos)
        Output: (gap_start_ind, gap_end_ind), tuple of 2 arrays:
                gap_start_ind - starting positions of the gaps
                gap_end_ind - ending positions of the gaps (including in the gap)
                arrays are arranged by increasing index of the gaps from the left end of the ref seq
    '''
    if len(gaps) == 0:
        return array([],dtype=int32), array([],dtype=int32)
    else:
        adjac = gaps[1:] - gaps[:-1] # difference between a position in the gap vec and the previous position, if the difference is 1, then it's consecutive gap
        gap_starts = where(adjac != 1)[0] # positions where new gap starts
        if len(gap_starts) == 0:
            gap_start_ind = array([gaps[0]])
            gap_end_ind = array([gaps[-1]])
        else:
            gap_start_ind = concatenate( ( array([gaps[0]]), gaps[ gap_starts + 1] )) # starting indices of all the consecutive gaps
            gap_end_ind = concatenate( ( gaps[ gap_starts ], array([gaps[-1]])))
        return gap_start_ind, gap_end_ind
## ======================================================================
def split_at_gap(gaps, ref_array):
    ''' Split the proposed PacBio sequence at the gap positions
        Input:  gaps - gap positions as an array, output from function gap_pos
                ref_array - 1d array for the proposed PacBio sequence
        Output: contiguous_seqs - list of contiguous parts of the PacBio sequence
    '''
    # return as a list of sequences, if there is no gap, the list will have length 1
    ref_short, ref_long = array_to_seq(ref_array)
    if len(gaps) == 0:
        return [ref_short]
    else:
        gap_start_ind, gap_end_ind = get_gaps(gaps)
        contiguous_seqs = []
        if gap_start_ind[0] != 0:
            seq = re.sub('-','',ref_long[:gap_start_ind[0]])
            contiguous_seqs.append(seq)
        for i in xrange(len(gap_start_ind) - 1):
            seq = re.sub('-','',ref_long[(gap_end_ind[i] + 1):gap_start_ind[i+1]])
            contiguous_seqs.append(seq)
        if gap_end_ind[-1] != (len(ref_long) - 1):
            seq = re.sub('-','',ref_long[(gap_end_ind[-1] + 1): len(ref_long)])
            contiguous_seqs.append(seq)
        return contiguous_seqs
## ======================================================================
def longest_seg(ref_array, read_array, minOverlap=10, minOverlapRatio=0.1):
    ''' Given a proposed PacBio sequence, find the length of the longest contiguous segment resulted from this ref array. The ref array might be so bad that it doesn't have any read support...
        Input:  ref_array, read_array - 1d and 2d array for the PacBio and all reads information
                minOverlap - minimum overlap length threshold
        Output: maxLen - length of the longest contiguous segment
    '''
    Cvec = get_compatible_reads(ref_array, read_array)
    if len(Cvec) == 0: # no read support at all
        return 0
    else:
        gaps = gap_pos(ref_array, read_array, Cvec, minOverlap=minOverlap,minOverlapRatio=minOverlapRatio)
        contiguous_seqs = split_at_gap(gaps, ref_array)
        contiguous_lengths = map(len, contiguous_seqs) # get length for each contiguous segment
        if len(contiguous_lengths) != 0:
            return max(contiguous_lengths)
        else:
            return 0
## ======================================================================
def get_reads_for_gap(read_array, gap, skip_reads=[]):
    ''' Given a gap (start, end), find reads that cover any base in the gap region and number of bases covered by them.
        Input:  read_array - array of read information
                gap - (start_pos, end_pos) of the gap trying to cover
                skip_reads - indices of reads to exclude (the ones that were selected from previous round for example)
        Output: reads_ind - indices of reads that cover at least 1 bp in the gap region
                reads_cov - corresponding numbers of bps covered by all the reads in reads_ind
    '''
    if len(skip_reads) > 0:
        keep_reads = array( [ i for i in arange(read_array.shape[0]) if i not in skip_reads ], dtype=int32 ) # indices of the reads that are kept (not skipped)
        read_array = read_array[keep_reads,:] # slicing the read_array
    else:
        keep_reads = arange(read_array.shape[0])

    if len(keep_reads) == 0:
        return array([],dtype=int32), array([],dtype=int32)
    else:
        gap_start = gap[0]
        gap_end = gap[1] + 1
        read_array = read_array[:, gap_start*5:gap_end*5 ] # only look at the particular positions
        coverages = apply_along_axis(cov_bps, axis=1, arr=read_array) # apply function cov_bps to find the coverage of all the other reads
        reads_ind = keep_reads[ where(coverages!=0)[0] ] # indices of reads that have nonzero coverage in this region
        reads_cov = coverages[ where(coverages!=0)[0] ] # their corresponding coverage
        return reads_ind, reads_cov # return both the indices and the corresponding coverage of the specified region
## ======================================================================
def get_new_ref(ref_array, read_ind, read_array):
    ''' from the previous ref_array, and the newly added read_array1d to fill the gap, find a new ref_array
        Input:  read_ind - read based on which the new ref_array is determined
                ref_array, read_array - as usual
        Output: ref1 - new ref_array with the newly added read incorporated
    '''
    ref1 = ref_array.copy()
    r_array= read_array[read_ind].reshape(-1,5) # convert the read's information to a 2d array with 5 columns
    cov_pos = unique( where( r_array == 1) [0] ) # find the positions covered by this read
    #print cov_pos
    # check each position: TODO now it's done with a loop, maybe should modify to a better way later TODO
    for pos in cov_pos: #TODO: actually, is it better to check only the region of the gap to be filled?
        read = r_array[pos,] # for each position covered by the new read, see if it's compatible with old ref
        ref = ref1[pos*5:(pos+1)*5]
        if not is_compatible(read, ref): # if the length 5 vectors are not compatible with each other
            #sys.stdout.write("At pos {} ref call {} changes to".format(pos, str(ref))) # print message about base changing in the new ref
            ref1[pos*5:(pos+1)*5] = 0
            read_calls = where(read == 1)[0]
            #print read_calls
            if len(read_calls) == 1: # only 1 call
                ref1[pos*5 + read_calls] = 1
            else: # if there are more than 1 option, pick the one with more read support
                #print "more than 1 options!!"
                pos_sum = read_array[ :, (pos*5 + read_calls)].sum(axis=0) # find the column sum for those calls
                pos_max = argmax(pos_sum)
                ref1[pos*5 + pos_max] = 1
            # sys.stdout.write(" {}\n".format(str(ref1[pos*5:(pos+1)*5]))) # for DEBUG
    return ref1
## ======================================================================
def write_compatible_reads(readsFasta, readinfo, compatible_ind, outDir):
    ''' Write read sequences in a fasta file.
        Input:  readsFasta - dictionary with read sequence information read_name => read_sequence
                readinfo - dictionary read_string => names of reads that have the read string
                compatible_ind - indices of reads that are compatible with some ref sequence
                outDir - output directory
        Output: fasta file in the output directory with read sequences
    '''
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    reads_names = get_reads_name(readinfo, compatible_ind)
    if len(reads_names) > 0: # only if there are reads provided to be printed
        reads_out = open(outDir + "/reads.fasta", 'w')
        for name in reads_names:
            if name in readsFasta:
                reads_out.write('>{}\n{}\n'.format(name, readsFasta[name]))
        reads_out.close()
    else:
        sys.stdout.write("compatible read list is empty, nothing to write.\n")
## ======================================================================
def greedy_fill_gap_1(read_array, ref0=None, minOverlap=10,minOverlapRatio=0.1, verbose=False):
    ''' Try to fill THE widest gap(just one gap, not all gaps) resulted from ref0 and minimize the number of uncovered bases using greedy algorithm
        If in verbose mode, write the newly proposed sequence and its compatible reads to files in a directory
        Input:  read_array - array including read information
                ref0 - ref_array to start with, if not specified, call the consensus sequence instead
        Output: (ref1, tot_gap, Cvec) - (new improved ref1, corresponding sequence with nucleotides, total number of gap positions/length, compatible read index array)
                files that include new PacBio sequence and its compatible reads (fasta)
    '''
    # if no sequence to start with, only call the consensus sequence, do not try to fill the gap
    if ref0 is None:
        if verbose:
            sys.stdout.write("Initial sequence not specified, use the read array to find consensus sequence and its gaps instead. \n")
        ref0 = get_consensus_from_array(read_array) # start with consensus sequence, summarized from all the reads
        Cvec = get_compatible_reads(ref0, read_array) # indices of reads that are compatible with ref0
        Gap_pos = gap_pos(ref0, read_array, Cvec, minOverlap, minOverlapRatio=minOverlapRatio) # positions not covered by the compatible reads (gap positions)

        ## starting sequence in string format for the consensus sequence
        seq0 = array_to_seq(ref0)[-1]
        for i in Gap_pos: # gap position will be written in lower case instead of upper case
            seq0 = seq0[:i] + seq0[i].lower() + seq0[(i+1):] 
        if verbose:
            sys.stdout.write("Consensus sequence is:\n{}\n".format(re.sub('-','',seq0)))

        if len(Gap_pos) > 0:
            gap_start_ind, gap_end_ind = get_gaps(Gap_pos) # starting and ending positions of all the gaps, in left to right order
            sys.stdout.write("\nGaps so far:\n")
            for i in xrange(len(gap_start_ind)):
                sys.stdout.write("({},{})\t".format(gap_start_ind[i], gap_end_ind[i]))
            sys.stdout.write("\n")
            gap_lens = gap_end_ind - gap_start_ind + 1 # gap lengths
            tot_gap = sum(gap_lens) # total gap size
            Mgap_ind = argmax(gap_lens) # index of the maximum gap among all gaps
            Mgap_len = gap_lens[Mgap_ind] # width of the maximum gap
            if verbose:
                sys.stdout.write("\n=== Maximum gap length is {}: ({}, {}).\n".format(Mgap_len, gap_start_ind[Mgap_ind], gap_end_ind[Mgap_ind]))
        else:
            tot_gap = 0
            if verbose:
                sys.stdout.write("\n=== No gap.\n")

        return ref0, tot_gap, Cvec
        #sys.stdout.write("   reads_ind: {}\n   reads_cov: {} \n".format(str(reads_ind), str(reads_cov))) # DEBUG

    # if a start sequence is provided, try to fill the gap introduced by this sequence with the following iterations
    # First try to fill the widest gap there, if it cannot be filled unless by increasing more gaps, move on to the next widest gap.
    # Stop whenever the total length of gaps decreases, or none of the gap could be filled at all.
    else:
        if verbose:
            sys.stdout.write("Initial sequence specified, now try to fill the gap in this sequence. \n")
        Cvec = get_compatible_reads(ref0, read_array) # indices of reads that are compatible with ref0
        Gap_pos = gap_pos(ref0, read_array, Cvec, minOverlap, minOverlapRatio) # positions not covered by the compatible reads (gap positions)

        # starting sequence in string format for the
        seq0 = array_to_seq(ref0)[-1] # with - for positions to delete
        for i in Gap_pos: # gap position will be written in lower case instead of upper case
            seq0 = seq0[:i] + seq0[i].lower() + seq0[i+1 :] 

        gap_start_ind, gap_end_ind = get_gaps(Gap_pos) # starting and ending positions of all the gaps, in left to right order
        gap_lens = gap_end_ind - gap_start_ind + 1 # gap lengths
        # sort gaps by their lengths
        ind_gap_sort = argsort(gap_lens)[::-1] # sort in increasing order and then reverse the array
        gap_start_ind = gap_start_ind[ind_gap_sort]
        gap_end_ind = gap_end_ind[ind_gap_sort]
        gap_lens = gap_lens[ind_gap_sort]
        
        if verbose:
            sys.stdout.write("\nGaps so far:\n")
            for i in xrange(len(gap_lens)):
                sys.stdout.write("({}, {})\t".format(gap_start_ind[i], gap_end_ind[i]))

        Min_tot_gap = sum(gap_lens) # smallest total gap size, initialize to the current total gap size
        Mgap_ind = ind_gap_sort[0] # index of the maximum gap among all gaps
        Mgap_len = gap_lens[0] # width of the maximum gap
        if verbose:
            sys.stdout.write("\n=== Maximum gap length is {}: ({}, {}).\n".format(Mgap_len, gap_start_ind[0], gap_end_ind[0]))

        best_ref = array(ref0)
        best_gap_pos = array(Gap_pos)
        best_Cvec = array(Cvec)
        improved = False
        gap_ind = 0

        while not improved and gap_ind <= (len(gap_lens) - 1): # until gap length improved, or all the gaps have been investigated
            if verbose:
                sys.stdout.write("gap_ind:  {} : ({}, {})\n".format(gap_ind,gap_start_ind[gap_ind], gap_end_ind[gap_ind]))
            reads_ind, reads_cov = get_reads_for_gap(read_array, (gap_start_ind[gap_ind], gap_end_ind[gap_ind]), skip_reads=Cvec) # get reads that can fill at least 1 base of the gap, and how many bases they fill
            # initialize the current best choice
            totally_filled = False # whether a gap is totally filled by current step, no gap filling for now, just initialization
            # start the iteration while loop
            # Stop condition:
            # 1. No more gap in the ref seq
            # 2. A gap is totally filled (ready to move on to the next gap)
            # 3. No more gap filling reads to try
            while len(best_gap_pos) > 0 and (not totally_filled) and len(reads_ind) > 0 :
                if verbose:
                    sys.stdout.write("gap positions: {}, remaining reads: {}, continue\n".format(len(best_gap_pos), len(reads_ind)))
                # sort the reads by the length of the gaps they fill
                ind_sort = argsort(reads_cov)[::-1] 
                reads_ind = reads_ind[ind_sort]
                reads_cov = reads_cov[ind_sort]
                ref1 = get_new_ref(ref0, reads_ind[0], read_array) # get a new ref, according to the highest ranked read

                Cvec1 = get_compatible_reads(ref1, read_array) # indices of reads that are compatible with ref1
                gap_pos1 = gap_pos(ref1, read_array, Cvec1, minOverlap, minOverlapRatio) # positions not covered by the compatible reads (gap positions)
                if len(gap_pos1) == 0:
                    if verbose:
                        sys.stdout.write("no more gaps\n")
                    Min_tot_gap = 0
                    best_ref = ref1
                    best_gap_pos = gap_pos1
                    improved = True
                    totally_filled = True
                    best_Cvec = Cvec1
                    break
                gap_start_ind1, gap_end_ind1 = get_gaps(gap_pos1) # starting and ending positions of all the gaps, in left to right order
                gap_lens1 = gap_end_ind1 - gap_start_ind1 + 1 # gap lengths
                if sum(gap_lens1) < Min_tot_gap:
                    improved = True
                    Min_tot_gap = sum(gap_lens1)
                    best_ref = ref1
                    best_gap_pos = gap_pos1
                    best_Cvec = Cvec1

                # if this step decreased the number of gaps by at least 1, and maximum gap is among them, then stop iteration
                if (len(gap_lens) - len(gap_lens1) == gap_lens[gap_ind]) and gap_start_ind[gap_ind] in setdiff1d(gap_start_ind, gap_start_ind1) and gap_end_ind[gap_ind] in setdiff1d(gap_end_ind, gap_end_ind1):
                    totally_filled = True
                    Min_tot_gap = 0
                    best_ref = ref1
                    best_gap_pos = gap_pos1
                    best_Cvec = Cvec1
                    improved = True

                remaining_inds = [] #  trying to find the indices in the array reads_ind that can still be tried to fill the gaps
                for r in reads_ind:
                    if r not in Cvec1: # delete the ones compatible with this chosen one, test these and see if the improvement is bigger. save remaining reads to check,
                        remaining_inds.append(where(reads_ind == r)[0][0])
                remaining_inds = array(remaining_inds)
                #print remaining_inds # DEBUG
                # update the list of remaining reads' information
                if len(remaining_inds) > 0 :
                    reads_ind = reads_ind[remaining_inds]
                    reads_cov = reads_cov[remaining_inds]
                else:
                    reads_ind = []
                    reads_cov = []
                #sys.stdout.write("\t   reads_ind: {}\n\t   reads_cov: {} \n\n".format(str(reads_ind), str(reads_cov))) # DEBUG
            if verbose:
                if len(best_gap_pos) == 0:
                    sys.stdout.write("no more gaps\n")
                elif totally_filled:
                    sys.stdout.write("gap totally filled\n")
                elif len(reads_ind) == 0:
                    sys.stdout.write("all reads tried\n")
                else:
                    sys.stdout.write("gap positions: {}, remaining reads: {}, continue\n".format(len(best_gap_pos), len(reads_ind)))
            gap_ind += 1

        if verbose:
            seq1 = array_to_seq(best_ref)[-1] # gap position will be written in lower case instead of upper case
            print best_gap_pos
            for i in best_gap_pos:
                seq1 = seq1[:i] + seq1[i].lower() + seq1[(i+1):] 
            sys.stdout.write("newly proposed sequence is:\n{}\n".format(re.sub('-','',seq1)))
            sys.stdout.write('\n')
            print_seqs(seq0, seq1) # show the difference between the starting sequence and the ending sequence after gap filling
    #sys.stdout.write("\n")
        return best_ref, Min_tot_gap, best_Cvec
## ======================================================================
def greedy_fill_gap(read_array, ref0=None, minOverlap=10, minOverlapRatio=0.1, verbose=False):
    ''' Try to fill THE widest gap(just one gap, not all gaps) resulted from ref0 and maximize the length of the longest segment in the resutled PacBio sequence using greedy algorithm
        If in verbose mode, write the newly proposed sequence and its compatible reads to files in a directory
        Input:  read_array - array including read information
                ref0 - ref_array to start with, if not specified, call the consensus sequence instead
        Output: (ref1, tot_gap, Cvec) - (new improved ref1, corresponding sequence with nucleotides, total number of gap positions/length, compatible read index array)
                files that include new PacBio sequence and its compatible reads (fasta)
    '''
    # if no sequence to start with, use consensus sequence as starting sequence instead
    if ref0 is None:
        if verbose:
            sys.stdout.write("Initial sequence not specified, use consensus sequence as starting sequence instead. \n")
        ref0 = get_consensus_from_array(read_array) # start with consensus sequence, summarized from all the reads
        
        #for i in xrange(read_array.shape[0]):
        #    print "read ", i+1
        #    print array_to_seq(read_array[i,:])[-1] # gap position will be written in lower case instead of upper case
        #    print "\n"

        #print "consensus sequence:"
        #print array_to_seq(ref0)[-1] # gap position will be written in lower case instead of upper case

    # if a start sequence is provided, try to fill the gap introduced by this sequence with the following iterations
    # Stop whenever the length of longest contiguous segment increases, or this length won't increase at all
    # First try to fill gaps from the widest to the smallest
    if verbose:
        sys.stdout.write("Try to fill the gap in this sequence. \n")
    Cvec = get_compatible_reads(ref0, read_array) # indices of reads that are compatible with ref0
    Gap_pos = gap_pos(ref0, read_array, Cvec, minOverlap, minOverlapRatio) # positions not covered by the compatible reads (gap positions)
    # starting sequence in string format
    seq0 = array_to_seq(ref0)[-1] # with - for positions to delete
    if len(Gap_pos) == 0: # If already covering the whole sequence
        return ref0, len(seq0)-seq0.count('-'), Gap_pos, Cvec

    ## If there is still room to improve the sequence
    for i in Gap_pos: # gap position will be written in lower case instead of upper case
        seq0 = seq0[:i] + seq0[i].lower() + seq0[(i+1):] 
    gap_start_ind, gap_end_ind = get_gaps(Gap_pos) # starting and ending positions of all the gaps, in left to right order
    gap_lens = gap_end_ind - gap_start_ind + 1 # gap lengths
    # sort gaps by their lengths
    ind_gap_sort = argsort(gap_lens)[::-1] # sort in increasing order and then reverse the array
    gap_start_ind = gap_start_ind[ind_gap_sort] # sort the gap indices and the lengths in decreasing order of the gap length
    gap_end_ind = gap_end_ind[ind_gap_sort]
    gap_lens = gap_lens[ind_gap_sort]
    
    # initialize the best solutions, and switches of ending conditions
    best_ref = ref0.copy()
    best_gap_pos = Gap_pos.copy()
    best_Cvec = Cvec.copy()
    try:
        old_max_len =  max(map(len, split_at_gap(Gap_pos, ref0)))
    except ValueError:
        old_max_len = 0

    best_max_len = old_max_len
    improved = False
    gap_ind = 0

    if verbose: # log message
        sys.stdout.write("\nGaps in starting sequence:\n")
        for i in xrange(len(gap_lens)):
            sys.stdout.write("({}, {}): {}\t".format(gap_start_ind[i], gap_end_ind[i], gap_lens[i]))
        sys.stdout.write("\nMaximum length is: {}\n".format(old_max_len))
        sys.stdout.write("=== Maximum gap length is {}: ({}, {}).\n".format(gap_lens[0], gap_start_ind[0], gap_end_ind[0]))


    while not improved and gap_ind < len(gap_lens): # until gap length improved, or all the gaps have been investigated
        if verbose:
            sys.stdout.write("Working on gap_ind:  {} : ({}, {})\n".format(gap_ind,gap_start_ind[gap_ind], gap_end_ind[gap_ind]))
        reads_ind, reads_cov = get_reads_for_gap(read_array, (gap_start_ind[gap_ind], gap_end_ind[gap_ind]), skip_reads=Cvec) # get reads that can fill at least 1 base of the gap, and how many bases they fill
        ind_sort = argsort(reads_cov)[::-1]  # sort the read indices by the coverage of the gap
        reads_ind = reads_ind[ind_sort]
        # initialize the current best choice
        totally_filled = False # whether a gap is totally filled by current step, no gap filling for now, just initialization
        # start the iteration while loop
        # Stop condition:
        # 1. No more gap in the ref seq
        # 2. A gap is totally filled (ready to move on to the next gap)
        # 3. No more gap filling reads to try
        while len(best_gap_pos) > 0 and (not totally_filled) and len(reads_ind) > 0 :
            # sort the reads by the length of the gaps they fill
            ref1 = get_new_ref(ref0, reads_ind[0], read_array) # get a new ref, according to the highest ranked read
            #print array_to_seq(ref0)[-1] # gap position will be written in lower case instead of upper case
            #print array_to_seq(ref1)[-1] # gap position will be written in lower case instead of upper case
            Cvec1 = get_compatible_reads(ref1, read_array) # indices of reads that are compatible with ref1
            gap_pos1 = gap_pos(ref1, read_array, Cvec1, minOverlap, minOverlapRatio) # positions not covered by the compatible reads (gap positions)
            try:
                max_len = max(map(len, split_at_gap(gap_pos1, ref1)))
            except ValueError:
                max_len = 0
            if len(gap_pos1) == 0: # by luck, all the gaps are filled!
                if verbose:
                    sys.stdout.write("no more gaps\n")
                improved = True
                best_max_len = max_len
                best_ref = ref1.copy()
                best_gap_pos = gap_pos1.copy()
                totally_filled = True
                best_Cvec = Cvec1.copy()
            else: # see if this one gap is totally filled
                gap_start_ind1, gap_end_ind1 = get_gaps(gap_pos1) # starting and ending positions of all the gaps, in left to right order
                gap_lens1 = gap_end_ind1 - gap_start_ind1 + 1 # gap lengths
                # if this step decreased the number of gaps by at least 1, and maximum gap is among them, then stop iteration
                if len(setdiff1d(gap_pos1,Gap_pos)) == 0 and (sum(gap_lens1) - sum(gap_lens)) == gap_lens[gap_ind]:
                    if verbose:
                        sys.stdout.write("This one gap is totally filled\n")
                    totally_filled = True
                    best_max_len = max_len
                    best_ref = ref1.copy()
                    best_gap_pos = gap_pos1.copy()
                    best_Cvec = Cvec1.copy()
                    improved = True
                else: # see if the gap is at least filled without making the contiguous read into pieces
                    if max_len > best_max_len:
                        improved = True
                        best_max_len = max_len
                        best_ref = ref1.copy()
                        best_gap_pos = gap_pos1.copy()
                        best_Cvec = Cvec1.copy()
                    reads_ind = delete(reads_ind,0)
                    delete_ind = []
                    for i in xrange(len(reads_ind)): # delete the ones compatible with this chosen one, test these and see if the improvement is bigger. save remaining reads to check,
                        if reads_ind[i] in Cvec1:
                            delete_ind.append(i)
                    reads_ind = delete(reads_ind,array(delete_ind))
                    if verbose:
                        sys.stdout.write("maximum length: {}, remaining reads: {}, continue\n".format(best_max_len, len(reads_ind)))
                #print array_to_seq(best_ref)[-1] # gap position will be written in lower case instead of upper case
        gap_ind += 1
        if verbose:
            sys.stdout.write("Next gap index: {}\n".format(gap_ind))

    if verbose:
        if improved:
            sys.stdout.write("\nGaps in resulted sequence:\n")
            best_start_ind, best_end_ind = get_gaps(best_gap_pos) # starting and ending positions of all the gaps, in left to right order
            best_lens = best_end_ind - best_start_ind + 1 # gap lengths
            for i in xrange(len(best_lens)):
                sys.stdout.write("({}, {}): {}\t".format(best_start_ind[i], best_end_ind[i], best_lens[i]))
            seq1 = array_to_seq(best_ref)[-1] # gap position will be written in lower case instead of upper case
            sys.stdout.write("\nImproved maximum contiguous length from {} to {}\n".format(old_max_len, best_max_len))
            for i in best_gap_pos:
                seq1 = seq1[:i] + seq1[i].lower() + seq1[(i+1):]
            sys.stdout.write('\n')
            print_seqs(seq0, seq1) # show the difference between the starting sequence and the ending sequence after gap filling
        else:
            sys.stdout.write("No improvement can be made\n")

#sys.stdout.write("\n")
    return best_ref, best_max_len, best_gap_pos, best_Cvec
## ======================================================================
def fill_gap(read_array, minOverlap=10, minOverlapRatio = 0.1, outFastaFile=None, outDir=None, readinfo=None, verbose=False):
    ''' Starting from the read_array, try to fill the gaps step by step until the result cannot be improved.
        Optionally, one can choose to output the PacBio sequence with smaller and smaller gaps in a file.

        Input:  read_array - array with read information
                minOverlap - minimum overlap length for reads to be considered
                outFastaFile - fasta file with all the reads originally mapped the PacBio sequence (needed for debug)
                outDir, readinfo - required for verbose mode
        Output: (best_ref, Max_len, Cvec) - (best ref_array so far, number of gap from this ref_array, compatible reads' indices)
    '''
    ref0 = greedy_fill_gap(read_array, ref0=None, minOverlap = minOverlap, minOverlapRatio=minOverlapRatio, verbose=verbose) # start with consensus sequence, summarized from all the reads
    iter_number = 1
    Max_len = ref0[1]
    print_tmp_files = False
    SeqLen = read_array.shape[1]/5

    if verbose:
        sys.stdout.write("verbose mode, write intermediate PacBio sequences and their compatible reads in files.\n")
        if outFastaFile is None:
            sys.stdout.write("Fasta file with mapped reads to Original PacBio sequence is not specified.\n")
        elif readinfo is None:
            sys.stdout.write("readinfo dictionary is not specified.\n")
        else:
            print_tmp_files = True
            if outDir is None:
                outDir = "./tmp/"
            if not os.path.exists(outDir): # make sure the output directory exists
                os.makedirs(outDir)
            if outDir[-1] != '/':
                outDir += '/'
            readsFasta = read_fasta(outFastaFile) # load read information to dictionary in memory

    while True:
        if print_tmp_files:
            sys.stdout.write("iteration number {}\n".format(iter_number))
            ref0_seq = array_to_seq(ref0[0])[0]
            cur_dir = outDir + '/round' + str(iter_number)
            if not os.path.exists(cur_dir):
                os.makedirs(cur_dir)
            seqOut = open(cur_dir + '/seq.fasta','w') # error corrected PacBio sequence at this point
            seqOut.write('>{} length: {} maximum length: {} compatible reads: {}\n{}\n'.format(iter_number, len(ref0_seq),  Max_len, len(ref0[-1]), ref0_seq))
            seqOut.close()
            samOut = open(cur_dir + '/compatible.sam','w') # corresponding sam file with compatible reads
            samfile_gen(read_array, ref0[0], readinfo, iter_number, len(ref0_seq), samOut, Cvec = ref0[-1])
            samOut.close()
            #write_compatible_reads(readsFasta, readinfo, ref0[-1], cur_dir + '/' )
        #print "ref0:", ref0
        if len(ref0[2]) > 0:
            ref1 = greedy_fill_gap(read_array, ref0[0], minOverlap, minOverlapRatio, verbose=verbose)
            max_len = ref1[1]
            if max_len > Max_len:
                Max_len = max_len
                ref0 = ref1
                iter_number += 1
                #compatible_reads = get_compatible_reads(ref0[0], read_array)
                #sys.stdout.write("number of compatible reads is {}, and the total number of reads is {}\n".format(len(compatible_reads), read_array.shape[0]))
            else:
                sys.stdout.write("not improving any more... Maximum contiguous length: {}, gap length: {}\n".format(Max_len, len(ref0[2])))
                break
        else:
            sys.stdout.write("no more gaps, finish iteration\n")
            break

    ## DEBUG -- figure out how many reads are compatible and non-compatible with the resulted PacBio sequence
    if verbose:
        compatible_reads = get_compatible_reads(ref0[0], read_array)
        sys.stdout.write("number of compatible reads is {}, and the total number of reads is {}\n".format(len(compatible_reads), read_array.shape[0]))
    return ref0
## ======================================================================
### Functions to generate sam files (or alignment in text) with updated PacBio sequence ###
## ======================================================================
def sam_record_gen(read_array1d, ref, qname, rname):
    ''' Generate sam record when a read is aligned to a reference sequence.
        Input:  read_array1d - read in 1d array format
                ref - reference/PacBio sequence in 1d array format
                qname - query/read name
                rname - reference sequence name
        Output: rec - alignment record in sam format
    '''
    ref_calls_pos = where(ref == 0)[0] # base pair positions where refrence sequence does not call
    read_array = array(read_array1d)
    read_array[ref_calls_pos] = 0
    read_short, read_long = array_to_seq(read_array1d,miss_char='-')
    ref_short, ref_long = array_to_seq(ref)
    ## generate CIGAR string and find leftmost mapping position on the long sequence
    CIGAR, rstart, rend = get_cigar(ref_long, read_long) 
    FLAG = '0' 
    if rstart == 0 :
        POS = 1
    else:
        POS = rstart - ref_long[:rstart].count('-') + 1
    POS = str(POS)
    MAPQ = '255'
    RNEXT = '='
    PNEXT = '0' 
    TLEN = '0'
    SEQ = re.sub('-','',read_short)
    QUAL = '*'
    char = re.findall('\D',CIGAR) # operation characters, MIDNSHP=X
    char = [x.upper() for x in char]  # convert to upper case
    count = map(int,re.findall('\d+',CIGAR))
    NM = str(sum(count[i] for i in xrange(len(char)) if char[i] in 'XID')) # length of the mapped read (including clipped part)
    rec = '\t'.join([qname, FLAG, str(rname), POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL, 'NM:i:'+NM]) + '\n'
    return rec
## ======================================================================
def samfile_gen(read_array, ref, readinfo, rname, rLen, out, Cvec=None):
    ''' Generate sam file given information of reads and reference in arrays.
        Input:  read_array - reads' information in array
                ref - reference/PacBio information in 1d array
                readinfo - dictionary read_string -> read_name
                rLen - length of the reference sequence
                rname - reference sequence name
                Cvec - indices of read_array's rows to include in sam file
        Output: sam file with the alignment of reads to the new error-corrected PacBio sequence
    '''
    if Cvec is None:
        Cvec = xrange(read_array.shape[0])
    #print Cvec
    read_names = [readinfo[key] for key in sorted(readinfo)]
    out.write('{}\t{}\t{}\n'.format('@HD', 'VN:1.4', 'SO:unsorted'))
    out.write('@SQ\tSN:{}\tLN:{}\n'.format(rname, rLen))
    out.write('@PG\tID:MetalRec\tPN:MetalRec\n')
    for i in Cvec:
        read_array1d = read_array[i,:]
        #print i, read_array1d
        qnames = read_names[i] # there could be more than one read with the same read info
        rec = sam_record_gen(read_array1d, ref, qnames[0], rname)
        out.write(rec)
        if len(qnames) >= 1:
            rec_list = rec.split('\t')
            for i in xrange(1, len(qnames)):
                rec_list[0] = qnames[i]
                out.write('\t'.join(rec_list))
