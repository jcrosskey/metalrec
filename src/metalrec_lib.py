#!/usr/bin/python

''' Collection of functions used in MetaLREC project
'''
import re
import sys, os
from numpy import *
import samread
from Bio import pairwise2 # pairwise alignment using dynamic programming
from Bio.pairwise2 import format_alignment
alphabet = 'ACGTD'
revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1]) # find reverse complement of a DNA sequence
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

        if reverse: # contigSeq => contigID
            # what if two proteins have the same sequences??
            # store the corresponding IDs in a list
            contigSeq = re.sub("\n","",contigSeq) # remove unix style EOF
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
    seq = re.sub('\n','',seq)
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
    try:
        flag = int(fields[1])
    except ValueError:
        sys.stderr.write('alignRecord is bad: \n {} \n'.format(alignRecord))
        return True

    if flag & 0x4 == 0x4: # read is unmapped (most reliable place to tell if a read is unmapped)
        return True
    elif cigarstring == '*': # if cigar string is not available, treat as bad 
        return True
    else: # read is NOT unmapped and cigar string is available
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

                    # For DEBUG
                    #if myread.qname == 'HISEQ11:283:H97Y1ADXX:1:1108:1971:87341':
                    #    #for i in realign_res:
                    #    #    print format_alignment(*i)
                    #    print "\n\n\n"
                    #    print format_alignment(*new_align)
                    #    print "\n\n\n"

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
def read_and_process_sam(samFile,rseq, maxSub=3, maxIns=3, maxDel=3,maxSubRate=0.02, maxInsRate=0.2, maxDelRate=0.2, minPacBioLen=1000, minCV=10,outsam='',outFastaFile=''):
    ''' Get consensus sequence from alignments of short reads to a long read, in the process, filter out bad reads and improve mapping

        Input:  samFile - sam file generated by mapping
                rseq - reference sequence as a string
                maxSub, maxIns, maxDel, maxSubRate, maxInsRate, maxDelRate are used to filter out badly mapped reads
                minPacBioLen - contiguous region length threshold to be a good region
                minCV - minimum coverage depth for a position to be considered as part of a good region
                outsam - output processed sam file name
                outFastaFile - fasta file including all the Illumina reads whose mappping passes the threshold

        Output: ref_bps - list of lists, one for each position on the reference sequence (without padding)
                ref_ins_dict - dictionary of insertions, one for each insertion position found in the alignments.
                readinfo - dictionary of read information, info_string => list of names of reads whose mapped segment corresponds to this info_string
    '''
    alphabet = 'ACGTD' # all the possible base pairs at a position
    keepRec = 0
    discardRec = 0
    lineNum = 0
    if outsam != '':
        newsam = open(outsam,'w')
    
    if outFastaFile != '':
        outFasta = open(outFastaFile, 'w')
    with open(samFile, 'r') as mysam:
        for line in mysam:
            lineNum += 1 
            #####################################
            # header line
            if line[0] == '@': 
                if outsam != '':
                    newsam.write(line) # if new reduced sam file is required, write the header lines into the new sam file
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
            #####################################
            # mapping record lines
            else:
                record = line.strip('\n')
                fields = record.split('\t')
                myread = samread.SamRead(record)
                if not is_record_bad(record, maxSub, maxIns, maxDel, maxSubRate, maxInsRate, maxDelRate): # if this alignment is good
                    keepRec += 1
                    if myread.is_paired():
                        print_name = myread.qname + ('/1' if myread.is_read1() else '/2')
                    else:
                        print_name = myread.qname

                    if outFastaFile != '':
                        if not myread.is_reverse():
                            outFasta.write('>{}\n{}\n'.format(print_name, myread.qSeq))
                        else:
                            outFasta.write('>{}\n{}\n'.format(print_name, revcompl(myread.qSeq)))

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

                    # if simplified sam file is required, find the new CIGAR string and write the new record
                    if outsam!= '':
                        cigarstring,first_non_gap = get_cigar(new_align[0], new_align[1]) # get the cigar string for the new alignment
                        # update information in the sam record
                        fields[3] = str(ref_region_start + new_align[3]) # starting position
                        fields[5] = cigarstring # cigar string
                        fields[9] = trimmed_qseq # trimmed query sequence
                        # write updated record in the new file
                        line = '\t'.join(fields) + '\n'
                        newsam.write(line)

                    # update string dictionary for the read information
                    read_string = dict_to_string(pos_dict) + ':' +  dict_to_string(ins_dict)
                    if read_string in readinfo:
                        readinfo[read_string].append(print_name)
                    else:
                        readinfo[read_string] = [print_name]

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
                    #if keepRec % 5000 == 0:
                    #    print lineNum
                    #    sys.stdout.write("discarded {} reads so far.\n".format(discardRec))
                    #    return ref_bps, ref_ins_dict
                else:
                    #print myread.qname
                    discardRec += 1
    
    if outsam != '':
        newsam.close()
    if outFastaFile != '':
        outFasta.close()
    sys.stdout.write("discarded {} reads.\n".format(discardRec))
    return ref_bps, ref_ins_dict, readinfo
## ======================================================================
def read_and_process_sam_pair(samFile,rseq, maxSub=3, maxIns=3, maxDel=3,maxSubRate=0.02, maxInsRate=0.2, maxDelRate=0.2, minPacBioLen=1000, minCV=10,outsam='',outFastaFile=''):
    ''' Get consensus sequence from alignments of short reads to a long read, in the process, filter out bad reads and improve mapping

        Input:  samFile - sam file generated by mapping
                rseq - reference sequence as a string
                maxSub, maxIns, maxDel, maxSubRate, maxInsRate, maxDelRate are used to filter out badly mapped reads
                minPacBioLen - contiguous region length threshold to be a good region
                minCV - minimum coverage depth for a position to be considered as part of a good region
                outsam - output processed sam file name, optional (discard reads that fail to pass the threshold, shuffle and clean up alignments)
                outFastaFile - fasta file including all the Illumina reads whose mappping passes the threshold, optional

        Output: ref_bps - list of lists, one for each position on the reference sequence (without padding, original positions in the PacBio sequence)
                ref_ins_dict - dictionary of insertions, one for each insertion position found in the alignments.
                readinfo - dictionary of read information, info_string => list of names of reads whose mapped segment corresponds to this info_string
    '''
    alphabet = 'ACGTD' # all the possible base pairs at a position
    keepRec = 0 # number of records kept and discarded
    discardRec = 0
    pairs = 0 # number of properly mapped pairs
    singles = 0 # number of reads that are mapped properly but not their mates
    lineNum = 0

    # checking and creating output files if specified
    if outsam != '':
        odir = os.path.dirname(os.path.abspath(outsam)) # make sure the directory for the output sam file exists
        if not os.path.exists(odir): # if not, create the directory
            os.makedirs(odir)
        newsam = open(outsam,'w')
    if outFastaFile != '':
        odir = os.path.dirname(os.path.abspath(outFastaFile)) # make sure the directory for the output fasta file exists
        if not os.path.exists(odir): # if not, create the directory
            os.makedirs(odir)
        outFasta = open(outFastaFile, 'w')
    
    mysam =  open(samFile, 'r') # input sam file object
    line = mysam.readline()
    is_pair = False # indicates if two reads are a properly mapped pair
    while 1: # start while loop to read all lines in the file
        if not line: # if at the end of the file, break out of the loop
            sys.stdout.write("File read complete.\n")
            break
        else: # otherwise, process the information in the file
            lineNum += 1 
            #####################################
            # header line
            if line[0] == '@': 
                if outsam != '':
                    newsam.write(line) # if new reduced sam file is required, write the header lines into the new sam file
                if line[1:3] == 'SQ': # reference sequence dictionary
                    rname = line[(line.find('SN:') + len('SN:')) : line.find('\t',line.find('SN:'))] # reference sequence name
                    rLen = int(line[(line.find('LN:') + len('LN:')) : line.find('\t',line.find('LN:'))]) # reference sequence length
                    sys.stdout.write("reference sequence {} has length {}. \n".format(rname, rLen))
                    if rLen < minPacBioLen:
                        sys.stdout.write("reference length is smaller than threshold {}.\n".format(minPacBioLen))
                        return 0 
                    ref_bps = [ [0] * 5 for x in xrange(rLen) ]  # list of lists, one list corresponding to every position on the reference sequence 'ACGTD' counts
                    ref_ins_dict = dict() # global insertion dictionary for the reference sequence
                    readinfo = dict() # dictionary storing read information (base call for the mapped read)
                line = mysam.readline() # continue reading next line
            #####################################
            # mapping record lines
            # Make sure that the paired end reads have the same name, instead arbitrarily considering two adjacent records to be a pair
            else:
                read_string = ''
                # First read of the pair, assuming mate reads are recorded next to each other in the sam file
                record_r1 = line.strip('\n')
                name_r1 = record_r1.split('\t')[0]
                # Try to see if next read is the Second read of the pair
                next_line = mysam.readline()
                if next_line == '':
                    is_pair = False
                else:
                    lineNum += 1
                    record_r2 = next_line.strip('\n')
                    name_r2 = record_r2.split('\t')[0]
                    if name_r1 == name_r2: # if their names match, then they are mates of each other's
                        is_pair = True
                    else:
                        is_pair = False

                # if r1 is mapped well
                if not is_record_bad(record_r1, maxSub, maxIns, maxDel, maxSubRate, maxInsRate, maxDelRate):
                    keepRec += 1 # number of kept records increases by 1

                    myread_r1 = samread.SamRead(record_r1)
                    if myread_r1.is_paired():
                        print_name = myread_r1.qname + ('/1' if myread_r1.is_read1() else '/2') #TODO: this might not be necessary ..
                    else:
                        print_name = myread_r1.qname
                        
                    if outFastaFile != '':
                        if not myread_r1.is_reverse():
                            outFasta.write('>{}\n{}\n'.format(print_name, myread_r1.qSeq))
                        else:
                            outFasta.write('>{}\n{}\n'.format(print_name, revcompl(myread_r1.qSeq)))
                    # re-map the read to the PacBio read
                    ref_region_start = max( myread_r1.rstart - 5, 1)
                    ref_region_end = min(myread_r1.get_rend() + 5, rLen)
                    trimmed_qseq = myread_r1.get_trim_qseq()
                    realign_res = pairwise2.align.globalms(rseq[(ref_region_start-1):ref_region_end], trimmed_qseq, 0, -1, -0.9, -0.9, penalize_end_gaps=[True, False])
                    new_align = pick_align(realign_res) # pick the best mapping: indel positions are the leftmost collectively,(seqA, seqB, score, first_non_gap_pos, last_non_gap_pos)
                    pos_dict, ins_dict = get_bases_from_align(new_align, ref_region_start + new_align[3])
                    # DEBUG
                    #if len(ins_dict) > 0:
                    #    print print_name
                    #    print format_alignment(*new_align)

                    # update corresponding base call frequencies for the reference 
                    for pos in pos_dict: # all the matching/mismatching/deletion positions
                        ref_bps[pos][alphabet.find(pos_dict[pos])] += 1 
                    for ins in ins_dict: # all the insertion positions
                        ins_chars = ins_dict[ins]
                        if ins not in ref_ins_dict: # if this position has not appeared in the big insertion dictionary, initialize it
                            ref_ins_dict[ins] = [[0,0,0,0] for ii in xrange(len(ins_chars))] # length is 4 because there is no 'D' at insertion position
                        elif len(ins_chars) > len(ref_ins_dict[ins]): # if the inserted bps is longer than the list corresponding to this position, make it longer
                            ref_ins_dict[ins] += [[0,0,0,0] for ii in xrange(len(ref_ins_dict[ins]),len(ins_chars))]
                        for i in xrange(len(ins_chars)):
                            ref_ins_dict[ins][i][alphabet.find(ins_chars[i])] += 1

                    # write the new mapping to the new/reduced sam file
                    fields = record_r1.split('\t')
                    # if simplified sam file is required, find the new CIGAR string and write the new record
                    if outsam!= '':
                        cigarstring,first_non_gap = get_cigar(new_align[0], new_align[1]) # get the cigar string for the new alignment
                        # update information in the sam record
                        fields[3] = str(ref_region_start + new_align[3]) # starting position
                        fields[5] = cigarstring # cigar string
                        fields[9] = trimmed_qseq # trimmed query sequence
                        # write updated record in the new file
                        line = '\t'.join(fields) + '\n'
                        newsam.write(line)

                    # update the read_string
                    if read_string != '':
                        read_string = read_string.split(':')[0] + dict_to_string(pos_dict) + ':' + read_string.split(':')[1] + dict_to_string(ins_dict)
                    else:
                        read_string = dict_to_string(pos_dict) + ':' + dict_to_string(ins_dict)

                    if is_pair: # if the next read is paired with the current read
                        # now check if the mate is also mapped well on to the PacBio read
                        if not is_record_bad(record_r2, maxSub, maxIns, maxDel, maxSubRate, maxInsRate, maxDelRate): # if r2 passes the threshold too
                            pairs += 1 # both reads from the pair passed the threshold and kept in the recuded sam file
                            myread_r2 = samread.SamRead(record_r2)
                            #print pairs, myread_r1.qname, myread_r2.qname # verbose output
                            keepRec += 1 # number of kept records increases by 1
                            if myread_r2.is_paired():
                                print_name = myread_r2.qname + ('/1' if myread_r2.is_read1() else '/2')
                            else:
                                print_name = myread_r2.qname
                                
                            if outFastaFile != '':
                                if not myread_r2.is_reverse():
                                    outFasta.write('>{}\n{}\n'.format(print_name, myread_r2.qSeq))
                                else:
                                    outFasta.write('>{}\n{}\n'.format(print_name, revcompl(myread_r2.qSeq)))

                            ref_region_start = max( myread_r2.rstart - 5, 1)
                            ref_region_end = min(myread_r2.get_rend() + 5, rLen)
                            trimmed_qseq = myread_r2.get_trim_qseq()
                            realign_res = pairwise2.align.globalms(rseq[(ref_region_start-1):ref_region_end], trimmed_qseq, 0, -1, -0.9, -0.9, penalize_end_gaps=[True, False])
                            new_align = pick_align(realign_res) # pick the best mapping: indel positions are the leftmost collectively
                            pos_dict, ins_dict = get_bases_from_align(new_align, ref_region_start + new_align[3])

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

                            fields = record_r2.split('\t')
                            # if simplified sam file is required, find the new CIGAR string and write the new record
                            if outsam!= '':
                                cigarstring,first_non_gap = get_cigar(new_align[0], new_align[1]) # get the cigar string for the new alignment
                                fields[3] = str(ref_region_start + new_align[3]) # starting position
                                fields[5] = cigarstring # cigar string
                                fields[9] = trimmed_qseq # trimmed query sequence
                                line = '\t'.join(fields) + '\n'
                                newsam.write(line)
                            # it's possible that the two pair-end reads overlap in the mapping, read string need special care in this case
                            # The conflict will be taken care of later when building the array for the reads
                            if read_string != '': # non-insertion string and insertion string are separated by :
                                read_string = read_string.split(':')[0] + dict_to_string(pos_dict) + ':' + read_string.split(':')[1] + dict_to_string(ins_dict)
                            else:
                                read_string = dict_to_string(pos_dict) + ':' + dict_to_string(ins_dict)
                        else:
                            discardRec += 1
                    else:
                        singles += 1 # count of reads that does not belong to a properly mapped pair

                # if r1 didn't pass the good read threshold, process mate instead
                else:
                    #print myread_r1.qname
                    discardRec += 1
                    if is_pair:
                        if not is_record_bad(record_r2, maxSub, maxIns, maxDel, maxSubRate, maxInsRate, maxDelRate): # if r2 passes the threshold too
                            singles += 1
                            keepRec += 1 # number of kept records increases by 1
                            myread_r2 = samread.SamRead(record_r2)
                            # printed name of this read, with /1 appended if it's the first read of the pair, otherwise append /2
                            if myread_r2.is_paired():
                                print_name = myread_r2.qname + ('/1' if myread_r2.is_read1() else '/2')
                            else:
                                print_name = myread_r2.qname

                            if outFastaFile != '':
                                if not myread_r2.is_reverse():
                                    outFasta.write('>{}\n{}\n'.format(print_name, myread_r2.qSeq))
                                else:
                                    outFasta.write('>{}\n{}\n'.format(print_name, revcompl(myread_r2.qSeq)))

                            ref_region_start = max(myread_r2.rstart - 5, 1)
                            ref_region_end = min(myread_r2.get_rend() + 5, rLen)
                            trimmed_qseq = myread_r2.get_trim_qseq()
                            realign_res = pairwise2.align.globalms(rseq[(ref_region_start-1):ref_region_end], trimmed_qseq, 0, -1, -0.9, -0.9, penalize_end_gaps=[True, False])
                            new_align = pick_align(realign_res) # pick the best mapping: indel positions are the leftmost collectively
                            pos_dict, ins_dict = get_bases_from_align(new_align, ref_region_start + new_align[3])

                            fields = record_r2.split('\t')
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

                            # if simplified sam file is required, find the new CIGAR string and write the new record
                            if outsam!= '':
                                cigarstring,first_non_gap = get_cigar(new_align[0], new_align[1]) # get the cigar string for the new alignment
                                # update information in the sam record
                                fields[3] = str(ref_region_start + new_align[3]) # starting position
                                fields[5] = cigarstring # cigar string
                                fields[9] = trimmed_qseq # trimmed query sequence
                                # write updated record in the new file
                                line = '\t'.join(fields) + '\n'
                                newsam.write(line)

                            if read_string != '':
                                read_string = read_string.split(':')[0] + dict_to_string(pos_dict) + ':' + read_string.split(':')[1] + dict_to_string(ins_dict)
                            else:
                                read_string = dict_to_string(pos_dict) + ':' + dict_to_string(ins_dict)

                        else:
                            #print myread_r2.qname
                            discardRec += 1

                # update string dictionary for the read information
                if read_string != '' and  read_string in readinfo:
                    if is_pair:
                        readinfo[read_string].append(print_name[:-2]) # pair-end reads are both properly mapped, trim the name
                    else:
                        readinfo[read_string].append(print_name) # only one is properly mapped, do not trim the name
                elif read_string != '':
                    if is_pair:
                        readinfo[read_string] = [print_name[:-2]]
                    else:
                        readinfo[read_string] = [print_name]

                if is_pair: # if the two reads form a pair, continue to read the next line
                    line = mysam.readline()
                else: # else need to process the second read separately
                    line = next_line
            if keepRec > 0 and keepRec % 1000 == 0:
                sys.stdout.write('  processed {} good records\n'.format(keepRec))
    mysam.close()

    if outsam != '':
        newsam.close()
    if outFastaFile != '':
        outFasta.close()
    sys.stdout.write("discarded {} reads.\n".format(discardRec))
    sys.stdout.write("discovered {} pairs of paired-end reads.\n".format(pairs))
    sys.stdout.write("discovered {} non-paired reads.\n".format(singles))
    sys.stdout.write("total number of lines is {}\n".format(lineNum))
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
            if seqA[ ins_pos - 1 ] == base: # if the base to the left is the same, shift it
                l_base = re.search(base+'+$', seqA[:ins_pos]).start()
                seqA = seqA[:l_base] + '-' + seqA[(l_base+1):ins_pos] + base + seqA[(ins_pos+1):]
        else:
            base = seqA[ins_pos] # base call corresponding to the insertion position
            if seqB[ ins_pos - 1 ] == base: # if the base to the left is the same, shift it
                l_base = re.search(base+'+$', seqB[:ins_pos]).start()
                seqB = seqB[:l_base] + '-' + seqB[(l_base+1):ins_pos] + base + seqB[(ins_pos+1):]

    return (seqA, seqB, score, begin, end)
## ======================================================================
def pick_align(align_list):
    ''' From a list of equivalent alignments between 2 sequences using dynamic progamming (Needleman-Wunch), pick the one whose indel positions are the most left
        Input:  align_list - list of tuples output from Bio.pairwise2.globalXX
        Output: align - one of the tuple in the list (seqA, seqB, score, first_non_gap_pos, last_non_gap_pos), last two elements were begin and end originally
    '''
    leftmost_indel_pos = (1000000,1000000)
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
            bestalign = (seqA, seqB, score, 0, len(seqA))
            First_non_gap = first_non_gap
            Last_non_gap = last_non_gap
    bestalign = shift_to_left(bestalign) # shift again to make sure the indels in the homopolymers are at the leftmost, this is not the case when number of best alignments is bigger than MAX_ALIGNMENTS in the pairwise2 module
    #print format_alignment(*bestalign)
    return bestalign[0], bestalign[1], bestalign[2], First_non_gap, Last_non_gap # return the list of aligns whose indel positions are the leftmost
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
def get_good_regions(ref_bps, rSeq, minPacBioLen=1000, minCV=10):
    ''' From mapping results from Illumina reads to a PacBio read, find the regions of PacBio read that satisfy these conditions:
        1.  Every base pair in the region is covered by at least 10 reads
        2.  The contiguous region has length >= 1000 bps
        
        Input:  ref_bps - dictionary for the match/mismatch and the deletion positions
                rSeq - PacBio sequence
        Output: good_regions - list of tuples, each tuple contain the begin and end coordinates of a good region
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
def get_poly_pos(ref_bps, ref_ins_dict, region=None, minReads=3, minPercent=0.01):
    ''' Get the polymorphic positions where the insertion happened, in the specified region (>=1000bps and covered by >=10 reads)
        Input:  ref_bps, ref_ins_dict - output objects from read_and_process_sam, noninsertion and insertion information for all the positions from the Illumina reads
                region - (begin, end) tuple of a region to consider
        Output: poly_bps - polymorphic non-insertion positions
                poly_ins - polymorphic insertion positions
                consensus_bps - nonpolymorphic non-insertion positions, get the consensus
                consensus_ins - nonpolymorphic insertion positions, get the consensus
    '''
    if region is None: # region to consider
        region = (0, len(ref_bps))
    poly_bps = [] # list of tuples
    poly_ins = []
    consensus_bps = []
    consensus_ins = []
    alphabet = 'ACGTD' # all the possible base pairs at a position
    cvs = [] # list of coverage depths

    # consider all the positions in the specified region
    for pos in xrange(region[0], region[1]):
        # check the match/mismatch/deletion first
        cv = sum(ref_bps[pos]) # coverage depth
        cvs.append(cv)
        base_calls = [ alphabet[i] for i in xrange(5) if ref_bps[pos][i] >= minReads and ref_bps[pos][i] >= cv * minPercent ]
        if len(base_calls) > 1:
            poly_bps.append((pos, base_calls))
        elif len(base_calls) == 1:
            consensus_bps.append((pos, base_calls[0]))
        # TODO: what if non of the base satisfies the condition? should we use information from PacBio read?
        else: # no base has more than required number of coverage, what now?
            base_calls = [ alphabet[i] for i in xrange(5) if ref_bps[pos][i] > 0 ]
            # if they all agree (low CV region) then use the consensus base
            if len(base_calls) == 1:
                consensus_bps.append((pos, base_calls[0]))
            # if they do not agree, pick the consensus base
            else:
                consensus_bps.append((pos, alphabet[ref_bps[pos].index(max(ref_bps[pos]))]))

        # check the insertion
        if pos in ref_ins_dict: 
            for i in xrange(len(ref_ins_dict[pos])): # look at all the bases inserted at this position
                pos_list = ref_ins_dict[pos][i]
                base_calls = [ alphabet[j] for j in xrange(4) if pos_list[j] >= minReads and pos_list[j] >= cv * minPercent ]
                if len(base_calls) > 1:
                    poly_ins.append(([pos,i], base_calls))
                elif len(base_calls) == 1:
                    consensus_ins.append(([pos,i], base_calls[0]))
                #else do not consider there is insertion at this position

    return poly_bps, poly_ins, consensus_bps, consensus_ins, cvs
## ======================================================================
def ref_extension(poly_bps, poly_ins, consensus_bps, consensus_ins, rseq, region=None, print_width = 100):
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
    if region is None: # if region is not specified, get it from poly_bps and consensus_bps
        begin = min([i[0] for i in consensus_bps + poly_bps])
        end = max([i[0] for i in consensus_bps + poly_bps])
        region = (begin,end+1)
    else:
        begin = region[0]
        end = region[1]

    ins_pos = [ i[0][0] for i in poly_ins + consensus_ins ] # all the insertion positions, if there is more than 1 base inserted, the position will be repeated
    bp_pos = [i[0] for i in poly_bps + consensus_bps] # all the non-insertion positions

    # distinct insertion positions (if there is more than 1 bp inserted, the insert position will appear more than once in ins_pos
    ins_pos_uniq = list(set(ins_pos))
    ins_pos_uniq.sort()

    poly_bp_pos = [i[0] for i in poly_bps]
    consensus_bp_pos = [i[0] for i in consensus_bps]

    ins_pos_dict = dict()
    bp_pos_dict = dict()

    inserted_bases = 0
    current_pos = region[0]
    newSeq = '' # new sequence with consensus base call positions updated
    oldSeq = '' # old sequence with just dashes for insertions
    matching_string = '' # matching string to print between two sequences for better visualization
    for ins in ins_pos_uniq:
        for mypos in xrange(current_pos, ins):
            if mypos not in consensus_bp_pos: # position has polymorphism or doesn't have enough short read coverage
                newSeq += rseq[mypos]
                oldSeq += rseq[mypos]
                matching_string += '|'
            else:
                newSeq += consensus_bps[consensus_bp_pos.index(mypos)][1] # if there is consensus
                oldSeq += rseq[mypos]
                if newSeq[-1] == oldSeq[-1]:
                    matching_string += '|'
                else:
                    matching_string += ' '
            bp_pos_dict[mypos] = mypos + inserted_bases # old position => new position in the extended ref sequence
        # now look at the insertion position
        newSeq += '-' * ins_pos.count(ins)
        oldSeq += '-' * ins_pos.count(ins)
        matching_string += ' ' * ins_pos.count(ins)
        ins_pos_dict[ins] = ins + inserted_bases
        inserted_bases += ins_pos.count(ins) # update number of inserted bases
        current_pos = ins

    # positions after the last insertion
    for mypos in xrange(current_pos, region[1]):
        if mypos not in consensus_bp_pos: # position has polymorphism or doesn't have enough short read coverage
            newSeq += rseq[mypos]
            newSeq += rseq[mypos]
            matching_string += '|'
        else:
            newSeq += consensus_bps[consensus_bp_pos.index(mypos)][1] # if there is consensus
            oldSeq += rseq[mypos]
            if newSeq[-1] == oldSeq[-1]:
                matching_string += '|'
            else:
                matching_string += ' '
        bp_pos_dict[mypos] = mypos + inserted_bases # old position => new position in the extended ref sequence

    # print the alignment between the original sequence and the newly proposed sequence
    # for the two different sequences, print the coordinates on each row, for both of them.
    for i in xrange(len(newSeq) / print_width + 1):
        sys.stdout.write('{:>5}\t{}\t{:<5}\n'.format(begin + i*print_width, re.sub('D','-',newSeq[i*print_width:(i+1)*print_width]), min(begin+(i+1)*print_width,len(newSeq)) ))
        sys.stdout.write('{:>5}\t{}\t{:<5}\n'.format(' ', matching_string[i*print_width:(i+1)*print_width],' ' ))
        sys.stdout.write('{:>5}\t{}\t{:<5}\n\n'.format(begin + i*print_width - oldSeq[:i*print_width].count('-'), oldSeq[i*print_width:(i+1)*print_width], min(begin+(i+1)*print_width - oldSeq[:(i+1)*print_width].count('-'),len(oldSeq)-oldSeq.count('-'))))

    return newSeq, bp_pos_dict, ins_pos_dict
## ======================================================================
def update_pos_info(poly_bps, poly_ins, consensus_bps, consensus_ins, bp_pos_dict, ins_pos_dict):
    ''' Given the correspondence between the positions in the original sequence and the positions in the extended sequence, update the polymorphic and consensus information for both insertion and non-insertion positions.
        Input:  poly_bps, poly_ins, consensus_bps, consensus_ins, bp_pos_dict, ins_pos_dict
        Output: poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext
    '''
    poly_bps_ext = []
    poly_ins_ext = []
    consensus_bps_ext = []
    consensus_ins_ext = []
    for i in poly_bps:
        poly_bps_ext.append( (bp_pos_dict[i[0]], i[1]) )
    for i in poly_ins:
        poly_ins_ext.append( (ins_pos_dict[i[0][0]] + i[0][1] , i[1]) )
    for i in consensus_bps:
        consensus_bps_ext.append( (bp_pos_dict[i[0]], i[1]) )
    for i in consensus_ins:
        consensus_ins_ext.append( (ins_pos_dict[i[0][0]] + i[0][1] , i[1]) )
    return poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext
## ======================================================================
def make_type_array(poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext):
    ''' Construct an array that indicates the type each position. 0: consensus non-insertion; 1: consensus insertion; 2: polymorphic non-insertion; 3: polymorphic insertion.
        Input:  position information list (including poly_bps_ext, poly_ins_ext, consensus_bps_ext, consensus_ins_ext) for the extended sequence
        Output: an array with 0,1,2,3 as entries and with length equal to the good region
    '''
    poly_bps_pos = array( [i[0] for i in poly_bps_ext] )
    poly_ins_pos = array( [ i[0] for i in poly_ins_ext] )
    consensus_bps_pos = array( [i[0] for i in consensus_bps_ext] )
    consensus_ins_pos = array( [ i[0] for i in consensus_ins_ext] )
    type_array = zeros(max(hstack((poly_bps_pos, poly_ins_pos, consensus_bps_pos, consensus_ins_pos))) + 1, dtype=int32)
    if len(consensus_ins_pos) > 0:
        type_array[consensus_ins_pos] = 1
    if len(poly_bps_pos) > 0:
        type_array[poly_bps_pos] = 2
    if len(poly_ins_pos) > 0:
        type_array[poly_ins_pos] = 3
    return type_array
## ======================================================================
def make_ref_array(consensus_bps_ext, consensus_ins_ext, type_array):
    ''' Construct array for the reference sequence, with consensus positions' entries filled in.
        Input:  consensus information (consensus_bps_ext, consensus_ins_ext)
                type_array - array including type of each position
        Output: ref_array - array with size 5 times size of type_array, each position has 5 options (ACGTD)
    '''
    alphabet = 'ACGTD'
    ref_array = zeros( len(type_array) * 5, dtype = int32 )
    for i in where(type_array == 0)[0]:
        base = [ j[1] for j in consensus_bps_ext if j[0] == i][0]
        #print base,i
        #print ref_array
        ref_array[ i*5 + alphabet.index(base) ] = 1
    for i in where(type_array == 1)[0]:
        base = [ j[1] for j in consensus_ins_ext if j[0] == i][0]
        ref_array[ i*5 + alphabet.index(base) ] = 1
    return ref_array
## ======================================================================
def array_to_seq(seq_array):
    ''' Convert 0-1 1d array back to nucleotide sequence.
        Input:  seq_array - 1d array with 0 and 1 as entries. length has to be multiple of 5, there can only be one 1 every 5 positions from start
        Output: seq - corresponding DNA sequence
    '''
    if len(seq_array) % 5 != 0: # array length has to be a multiple of 5
        sys.exit("length of the array is not a multiple of 5!! \n")
    else:
        seqLen = len(seq_array) / 5 # sequence length
        seq_array = seq_array.reshape(-1,5) # reshape to an array with 5 columns
        seq = ["x"] * seqLen # initialize the sequence to return, as a list
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
        sys.stdout.write('{:>5}\t{}\t{:<5}\n\n'.format(' ', seq2[i*print_width:(i+1)*print_width], ' '))
## ======================================================================
def make_read_array1d(read_string, bp_pos_dict, ins_pos_dict, type_array, poly_bps, poly_ins, consensus_bps, consensus_ins):
    ''' Make 1d array for a particular read from its string (key of dictionary readinfo).
        Incorporate the base calling from other reads at the same position.

        Input:  read_string - read string (key of dictionary readinfo)
                bp_pos_dict, ins_pos_dict - correspondence between positions from original to extended ref sequence
                length - length of the reference sequence ( good region )
                poly_bps, poly_ins, consensus_bps, consensus_ins - base calling information for all the positions (non-extended sequence)

        Output: read_array1d - 1d array for this particular read string
    '''
    read_array1d = zeros( len(type_array) * 5, dtype=int32 ) # initialization
    bpstring = read_string.split(":")[0] # non-insertion information
    insstring = read_string.split(":")[1] # insertion position's information

    # all the positions of a read that lay in the specified region (covered by the dicts)
    map_positions = []
    for pos in map(int, re.findall('\d+',bpstring)):
        if pos in bp_pos_dict:
            map_positions.append(bp_pos_dict[pos])
    if len(map_positions) == 0:
        return read_array1d
    else:
        start_pos = min(map_positions)
        end_pos = max(map_positions)
        # get the base calling information for all the positions, with coordinates in the extended sequence
        poly_bps, poly_ins, consensus_bps, consensus_ins = update_pos_info(poly_bps, poly_ins, consensus_bps, consensus_ins, bp_pos_dict, ins_pos_dict)

        # If this read covers some non-insertion positions, do the following:
        if len(bpstring) > 0:
            positions = map(int, re.findall('\d+',bpstring)) # positions
            bases = re.findall('\D+',bpstring) # base calling
            for i in xrange(len(positions)):
                position = bp_pos_dict[positions[i]] # position in the extended sequence
                if type_array[position] == 0: # consensus non-inseriton position, read's call should be the same as the consensus call
                    c_bp = [l[1] for l in consensus_bps if l[0] == position][0]
                    read_array1d[ position*5 + alphabet.index(c_bp) ] = 1
                if type_array[position] == 2: # polymorphic non-insertion position, check if read's call is one of the possible calls
                    p_bps = [l[1] for l in poly_bps if l[0] == position][0]
                    if bases[i] in p_bps: # if read's call is one of the possible calls, do not change it
                        read_array1d[ position*5 + alphabet.index(bases[i]) ] = 1
                    else: # if read's call is not among the possible calls, it's considered as an error, so it could be either of the possible calls
                          # TODO: this is done in somewhat awkward way, try to fix this
                        for p_bp in p_bps:
                            read_array1d[ position*5 + alphabet.index(p_bp) ] = 1
                            
        # If this read covers some insertion positions, do the following:            
        positions = map(int, re.findall('\d+',insstring))
        bases = re.findall('\D+',insstring)
        for i in xrange(len(positions)): # for each base inserted at this position (there might be more than 1 base inserted)
            if positions[i] in ins_pos_dict: # if this position is decided to be an insertion position in the extended sequence, proceed. Otherwise, it's treated as an error
                position = ins_pos_dict[positions[i]]
                for j in xrange(len(bases[i])):
                    ins_position = position + j
                    if type_array[ins_position] == 1: # consensus insertion position
                        c_bp = [l[1] for l in consensus_ins if l[0] == ins_position][0]
                        read_array1d[ ins_position*5 + alphabet.index(c_bp) ] = 1

                    if type_array[ins_position] == 3: # polymorphic insertion position
                        p_bps = [l[1] for l in poly_ins if l[0] == ins_position][0]
                        if bases[i][j] in p_bps:
                            read_array1d[ ins_position*5 + alphabet.index(bases[i][j]) ] = 1
                        else:
                            for p_bp in p_bps:
                                read_array1d[ position*5 + alphabet.index(p_bp) ] = 1

        # check if this read missed any insertion position
        for ins_pos in hstack((where(type_array==3)[0], where(type_array==1)[0])):
            if ins_pos >= start_pos and ins_pos <= end_pos and (ins_pos not in positions): # if there is no insertion for this read at the insertion position, put a deletion
                read_array1d[ins_pos * 5 + 4 ] = 1

        return read_array1d
## ======================================================================
def make_read_array(readinfo, bp_pos_dict, ins_pos_dict, type_array, poly_bps, poly_ins, consensus_bps, consensus_ins):
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
        read_array[i,:] = make_read_array1d(read_string, bp_pos_dict, ins_pos_dict, type_array, poly_bps, poly_ins, consensus_bps, consensus_ins)
        read_counts[i] = len(readinfo[read_string])
        i += 1
    return read_array, read_counts
## ======================================================================
def is_compatible(array1, array2):
    ''' Find out if the base calls of 2 reads at a certain position are compatible or not, given the 0-1 vectors representing the base calls. For now, the length of the vectors is the same (5 for ACGTD).
        Input:  array1, array2 - length 5 vectors, each corresponds to the base call of a read at one position
        Output: boolean value - true if they are compatible, false if not
    '''
    same_call_pos = where(bitwise_and( array1, array2) == 1)[0]
    if len(same_call_pos) > 0:
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
    nread = read_array.shape[0]
    Cmat = zeros((nread, nread),dtype=int32)
    for i in xrange(nread):
        for j in xrange(i+1, nread):
            if are_reads_compatible(read_array[i],read_array[j]):
                Cmat[i,j] = 1
    return array(Cmat)
## ======================================================================
def cov_bps(read_array1d):
    ''' Find number of bases covered by an array '''
    r = read_array1d.reshape(-1,5)
    cov_pos = unique( where( r==1 )[0] ) # positions covered by this read
    return len(cov_pos)
## ======================================================================
def is_read_compatible(ref_array, read_array1d):
    ''' Given a reference array and an array for a read, determine if they are compatible or not.
        Input:  ref_array - vector of 0-1 for base calls on the reference sequence
                read_array1d - 1d array from one read
        Output: True if read is compatible, False otherwise
    '''
    r = read_array1d.reshape(-1,5)
    cov_pos = unique( where( r==1 )[0] ) # positions covered by this read
    if dot(ref_array, read_array1d) == len(cov_pos):
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
        if is_read_compatible(ref_array, read_array[i,:]):
            compatible_ind.append(i)
    return array(compatible_ind)
## ======================================================================
def get_reads_name(readinfo, compatible_ind):
    ''' Find the names of reads that are compatible with a given PacBio sequence, given the readinfo dictionary and the compatible array row indices.
        Input:  readinfo - dictionary with read information read_string => list of read names that have this read_string
                compatible_ind - row indices of the array compatible with the PacBio sequence
        Output: reads_name_list - list of names of reads corresponding to the compatible_ind
    '''
    reads_name_list = [] # initialize an empty list
    info_keys = sorted(readinfo.keys()) # sort the keys of readinfo
    for ind in compatible_ind:
        reads_name_list += readinfo[ info_keys[ind] ]
    return reads_name_list
## ======================================================================
def gap_pos(ref_array, read_array, compatible_ind):
    ''' Given a reference array and indices (of array) of compatible reads, find the gap positions, i.e. positions that are not covered by any read
        Input:  ref_array - vector of 0-1 for base calls on the reference sequence
                read_array - 2d array from all the reads
                compatible_ind - indices of compatible reads, corresponding to the read_array
        Output: gap_vec - vector with positions that have no coverage from the reads
    '''
    read_array = read_array[compatible_ind,:]
    base_pos = where( ref_array == 1)[0] # coordinates corresponding to the bases called by the reference
    base_cov = read_array[ : , base_pos].reshape(-1,len(base_pos))
    base_cov = base_cov.sum(axis=0) # pick the compatible rows and the correct columns, sum over the columns
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
def get_reads_for_gap(read_array, gap, skip_reads=[]):
    ''' Given a gap (start, end), find reads that cover any base in the gap region and number of bases covered by them.
        Input:  read_array - array of read information
                gap - (start_pos, end_pos) of the gap trying to cover
                skip_reads - indices of reads to exclude (the ones that were selected from previous round for example)
        Output: reads_ind - indices of reads that cover at least 1 bp in the gap region
                reads_cov - corresponding numbers of bps covered by all the reads in reads_ind
    '''
    if len(skip_reads) > 0:
        keep_reads = array( [ i for i in arange(read_array.shape[0]) if i not in skip_reads ] ) # indices of the reads that are kept (not skipped)
        read_array = read_array[keep_reads,:] # slicing the read_array
    else:
        keep_reads = arange(read_array.shape[0])

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
        Output:
    '''
    ref1 = ref_array.copy()
    r_array= read_array[read_ind].reshape(-1,5) # convert the read's information to a 2d array with 5 columns
    #print r_array
    cov_pos = unique( where( r_array == 1) [0] ) # find the positions covered by this read
    #print cov_pos
    # check each position: TODO now it's done with a loop, maybe should modify to a better way later TODO
    for pos in cov_pos:
        read = r_array[pos,]
        ref = ref1[pos*5:(pos+1)*5]
        if not are_reads_compatible(read, ref): # if the length 5 vectors are not compatible with each other
            #sys.stdout.write("At pos {} ref call {} changes to".format(pos, str(ref))) # print message about base changing in the new ref
            ref1[pos*5:(pos+1)*5] = 0
            read_calls = where(read == 1)[0]
            #print read_calls
            if len(read_calls) == 1: # only 1 call
                ref1[pos*5 + read_calls] = 1
            else: # if there are more than 1 option, pick the one with more read support
                print "more than 1 options!!"
                pos_sum = read_array[ :, (pos*5 + read_calls)].sum(axis=0) # find the column sum for those calls
                pos_max = argmax(pos_sum)
                ref1[pos*5 + pos_max] = 1
            # sys.stdout.write(" {}\n".format(str(ref1[pos*5:(pos+1)*5]))) # for DEBUG
    return ref1
## ======================================================================
def write_compatible_reads(readsFasta, readinfo, compatible_ind, outDir):
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    reads_names = get_reads_name(readinfo, compatible_ind)
    reads_out = open(outDir + "/reads.fasta", 'w')
    for name in reads_names:
        if '/' in name: # only one read of the pair was properly mapped
            reads_out.write('>{}\n{}\n'.format(name, readsFasta[name]))
        else: # both reads were properly mapped
            name1 = name + '/1'
            name2 = name + '/2'
            if name1 in readsFasta:
                reads_out.write('>{}\n{}\n'.format(name1, readsFasta[name1]))
            if name2 in readsFasta:
                reads_out.write('>{}\n{}\n'.format(name2, readsFasta[name2]))

    reads_out.close()

## ======================================================================
def greedy_fill_gap(read_array, ref0=None, readsFasta=None, readinfo=None):
    ''' Try to fill THE widest gap(just one gap, not all gaps) resulted from ref0 and minimize the number of uncovered bases using greedy algorithm
        If in DEBUG mode, write the newly proposed sequence and its compatible reads to files in a directory
        Input:  read_array - array including read information
                ref0 - ref_array to start with, if not specified, call the consensus sequence instead
                readsFasta - dictionary between read name and its sequence
                readinfo - dictionary between read_string and its corresponding read names
        Output: (ref1, Min_gap) - (new improved ref1, total number of gap positions/length)
                files that include new PacBio sequence and its compatible reads (fasta)
    '''
    # if no sequence to start with, only call the consensus sequence, do not try to fill the gap
    if ref0 is None:
        sys.stdout.write("Initial sequence not specified, use the read array to find consensus sequence and its gaps instead. \n")
        ref0 = get_consensus_from_array(read_array) # start with consensus sequence, summarized from all the reads

        Cvec = get_compatible_reads(ref0, read_array) # indices of reads that are compatible with ref0

        Gap_pos = gap_pos(ref0, read_array, Cvec) # positions not covered by the compatible reads (gap positions)

        # starting sequence in string format for the
        seq0 = array_to_seq(ref0)[-1]
        for i in Gap_pos:
            seq0 = seq0[:i] + seq0[i].lower() + seq0[i+1 :] 

        gap_start_ind, gap_end_ind = get_gaps(Gap_pos) # starting and ending positions of all the gaps, in left to right order
        gap_lens = gap_end_ind - gap_start_ind + 1 # gap lengths
        Min_gap = sum(gap_lens) # smallest gap size, initialize to the current gap size
        Mgap_ind = argmax(gap_lens) # index of the maximum gap among all gaps
        Mgap_len = gap_lens[Mgap_ind] # width of the maximum gap

        sys.stdout.write("\n=== Maximum gap length is {}: ({}, {}).\n".format(Mgap_len, gap_start_ind[Mgap_ind], gap_end_ind[Mgap_ind]))
        return ref0, Min_gap, Cvec
        #sys.stdout.write("   reads_ind: {}\n   reads_cov: {} \n".format(str(reads_ind), str(reads_cov))) # DEBUG

    # stop condition:
    # 1. No more gap in the ref seq
    # 2. A gap is totally filled (ready to move on to the next gap)
    # 3. No more gap filling reads to try
    else:
        ref0 = ref0[0]
        Cvec = get_compatible_reads(ref0, read_array) # indices of reads that are compatible with ref0
        Gap_pos = gap_pos(ref0, read_array, Cvec) # positions not covered by the compatible reads (gap positions)

        # starting sequence in string format for the
        seq0 = array_to_seq(ref0)[-1]
        for i in Gap_pos:
            seq0 = seq0[:i] + seq0[i].lower() + seq0[i+1 :] 

        gap_start_ind, gap_end_ind = get_gaps(Gap_pos) # starting and ending positions of all the gaps, in left to right order
        gap_lens = gap_end_ind - gap_start_ind + 1 # gap lengths
        Min_gap = sum(gap_lens) # smallest gap size, initialize to the current gap size
        Mgap_ind = argmax(gap_lens) # index of the maximum gap among all gaps
        Mgap_len = gap_lens[Mgap_ind] # width of the maximum gap

        reads_ind, reads_cov = get_reads_for_gap(read_array, (gap_start_ind[Mgap_ind], gap_end_ind[Mgap_ind]), skip_reads=Cvec) # get reads that can fill at least 1 base of the gap, and how many bases they fill
        sys.stdout.write("\n=== Maximum gap length is {}: ({}, {}).".format(Mgap_len, gap_start_ind[Mgap_ind], gap_end_ind[Mgap_ind]))
        totally_filled = False # whether a gap is totally filled by current step, no gap filling for now, just initialization
        best_ref = ref0
        best_gap_pos = Gap_pos
        Cvec1 = Cvec
        while len(Gap_pos) > 0 and (not totally_filled) and len(reads_ind) > 0 :
            # sort the read indices and the fill lengths in descending order of the fill length
            ind_sort = argsort(reads_cov)[::-1] 
            reads_ind = reads_ind[ind_sort]
            reads_cov = reads_cov[ind_sort]
            #sys.stdout.write("   reads_ind: {}\n   reads_cov: {} \n\n".format(str(reads_ind), str(reads_cov))) # DEBUG
            #print len(reads_ind)
            # get a new ref, according to the highest ranked read
            ref1 = get_new_ref(ref0, reads_ind[0], read_array) 

            Cvec1 = get_compatible_reads(ref1, read_array) # indices of reads that are compatible with ref1
            gap_pos1 = gap_pos(ref1, read_array, Cvec1) # positions not covered by the compatible reads (gap positions)
            if len(gap_pos1) == 0:
                Min_gap = 0
                best_ref = ref1
                break
            gap_start_ind1, gap_end_ind1 = get_gaps(gap_pos1) # starting and ending positions of all the gaps, in left to right order
            gap_lens1 = gap_end_ind1 - gap_start_ind1 + 1 # gap lengths
            if sum(gap_lens1) < Min_gap:
                Min_gap = sum(gap_lens1) # update current best gap size and the corresponding ref_array
                best_ref = ref1
                best_gap_pos = gap_pos1

            # if this step decreased the number of gaps by at least 1, and maximum gap is among them, then stop iteration
            if (len(gap_lens) - len(gap_lens1) >= 1) and gap_start_ind[Mgap_ind] in setdiff1d(gap_start_ind, gap_start_ind1):
                totally_filled = True

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

        seq1 = array_to_seq(best_ref)[-1]
        for i in best_gap_pos:
            seq1 = seq1[:i] + seq1[i].lower() + seq1[i+1 :] 
        sys.stdout.write('\n')
        print_seqs(seq0, seq1) # show the difference between the starting sequence and the ending sequence after gap filling
    #sys.stdout.write("\n")
    return best_ref, Min_gap, Cvec1
## ======================================================================
def fill_gap(read_array, outFastaFile=None, outDir=None, readinfo=None):
    ''' Starting from the read_array, try to fill the gaps step by step until the result cannot be improved.
        Optionally, one can choose to output the PacBio sequence with smaller and smaller gaps in a file.

        Input:  read_array - array with read information
                outFastaFile - fasta file with all the reads originally mapped the PacBio sequence (needed for debug)
                #seq_file - file to store the PacBio sequence produced by the function each step
        Output: (best_ref, Mingap, Cvec) - (best ref_array so far, number of gap from this ref_array, compatible reads' indices)
                outDir - if specified, directory to save all the intermediate files
                #seq_file - if specified, step-by-step PacBio sequences
    '''
    ref0 = greedy_fill_gap(read_array) # start with consensus sequence, summarized from all the reads
    iter_number = 1
    Mingap = ref0[1]
    print_tmp_files = False

    if __debug__:
        sys.stdout.write("DEBUG mode, write intermediate PacBio sequences and their compatible reads in files.\n")
        if outFastaFile is None:
            sys.stdout.write("Fasta file with mapped reads to Original PacBio sequence is not specified.\n")
        elif readinfo is None:
            sys.stdout.write("readinfo dictionary is not specified.\n")
        else:
            print_tmp_files = True
            if outDir is None:
                outDir = "./tmp/"
            readsFasta = read_fasta(outFastaFile)

    while True:
        if print_tmp_files:
            ref0_seq = array_to_seq(ref0[0])[0]
            cur_dir = outDir + str(iter_number)
            if not os.path.exists(cur_dir):
                os.makedirs(cur_dir)
            seqOut = open(cur_dir + '/seq.fasta','w')
            seqOut.write('>{} gap length: {} compatible reads: {}\n{}\n'.format(iter_number,  Mingap, len(ref0[-1]), ref0_seq))
            seqOut.close()
            write_compatible_reads(readsFasta, readinfo, ref0[-1], cur_dir + '/' )
        #print "ref0:", ref0
        if Mingap > 0:
            ref1 = greedy_fill_gap(read_array, ref0)
            gaps = ref1[1]
            if gaps < Mingap:
                Mingap = gaps
                ref0 = ref1
                iter_number += 1
            else:
                print "not improving any more :("
                break
        else:
            break

    return ref0[0], ref0[1]
