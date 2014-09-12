#!/home/cjg/bin/python

'''
muscleAlign

muscleAlign is the second program in the pipeline for protein family clustering, 
align sequences in all clusters by MUSCLE

Created by JJ Chai on 09/20/2013.
Copyright (c) 2013 JJ Chai (ORNL). All rights reserved.

'''

# Import Python modules
import sys, warnings, os
import time
import glob
from datetime import datetime
from subprocess import Popen, PIPE, check_call, STDOUT

# Import local modules
sys.path.append("/chongle/jj/01_PfClust/scripts")
import jj_utils

## =================================================================
## argument parser
import argparse

## Version
version_str = "0.0.2"
'''
    1. multi-process
    2. next: check files in the directory, and add ability to use options other than default in MUSCLE
'''

parser = argparse.ArgumentParser(description="Multiple sequence alignments with MUSCLE",
    prog = 'muscleAlign', #program name
    prefix_chars='-', # prefix for options
    fromfile_prefix_chars='@', # if options are read from file, '@args.txt'
    conflict_handler='resolve', # for handling conflict options
    add_help=True, # include help in the options
    formatter_class=argparse.ArgumentDefaultsHelpFormatter # print default values for options in help message
    )
## version control
parser.add_argument("--version", action="version",version='%(prog)s {}'.format(version_str))

## --quiet and --verbose are mutually exclusive
group = parser.add_mutually_exclusive_group()
group.add_argument("-v", "--verbose", action="store_false",help="verbose mode, more output")
group.add_argument("-q", "--quiet", action="store_true",help="quiet mode, as little output as possible")


## input fasta files and/or directories that include fasta files to be aligned
parser.add_argument("-f","--fasta",nargs="+",help="fasta file with sequences to be aligned",dest='fastafiles')
parser.add_argument("-d","--dir",nargs='+',help="directory that include all the fasta files",dest='inputdirs')

## output directory
parser.add_argument("-o", "--out",nargs='?',help="directory for output files",const='./',default='./',dest='outputdir')

## multi-process
parser.add_argument("-p","--multi-processes",nargs='?',help="number of processes used for alignment",
    const=1,default=1,dest='num_proc',type=int)

## =================================================================
## get muscle command used in shell and the location of the output log file
def get_muscle_cmd_and_log(inputfile,outputdir='./'):

    ## full path of muscle: change if there is a new version or location has changed
    muscle_full = "/home/cjg/bin/muscle"    

    ## output directory, make sure it ends with "/"
    outputdir = jj_utils.set_directory_path(outputdir) 

    ## if input file is aa.fa, output is aa_align.fa, log is aa_align.log
    basename = os.path.basename(inputfile)
    name, suffix = os.path.splitext(basename)
    outputfile = outputdir + name + ".align_fa"  # output file name
    logfile = outputdir + name + ".align_log" # log file name

    ## muscle command, use default option for now
    muscle_cmd = "%s -in %s -out %s -log %s" %(muscle_full,inputfile,outputfile,logfile)

    return(muscle_cmd,logfile)

## =================================================================
## multi-processing
def muscleAlign(args):
    
    ## output directory from options, default is "./"
    outputdir = os.path.abspath(args.outputdir)
    ## if output dir does not exist, make this directory
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    ## print output directory information
    #sys.stderr.write("\noutput directory is %s" %outputdir)
    #sys.stderr.write("fastafiles are {}".format(args.fastafiles))
    
    ## create list to store all fasta files to be processed
    if args.fastafiles:
        fa_files = args.fastafiles
    else:
        fa_files = []

    ## if directories are provided, look for all the .fa and .fasta files in each directory
    if args.inputdirs:
        for dir in args.inputdirs:
            dir = os.path.abspath(dir)
            fa_files.extend(glob.glob(dir + "/*.fa")) # add .fa files to the list of files
            fa_files.extend(glob.glob(dir + "/*.fasta")) # add .fasta files to the list of files

    ## if there is no input file,exit
    if len(fa_files) is 0:
        jj_utils.die("no file to align! use -h to see help.")
    ## if there is only 1 file to align, set the number of processes to 1
    elif len(fa_files) is 1:
        number_proc = 1
        sys.stderr.write("Only 1 file to process, number of processes is set to: %3d\n" %number_proc)
        ## call Popen to run MUSCLE from shell                                                                          
        muscle_cmd,muscle_log = get_muscle_cmd_and_log(fa_files[0],outputdir) 
        process = Popen(muscle_cmd,shell=True, stdout=PIPE, stderr=STDOUT)                                              
        ## write log file
        muscle_log_file = open(muscle_log, 'wb')                                             
        cmd_output = process.communicate()[0].strip()                                                           
        muscle_log_file.write(cmd_output)                                                                       
        muscle_log_file.close()           
    ## if there are more than 1 file to process
    else: 
        ## number of alignments to be done
        fa_count = len(fa_files)
        number_proc = args.num_proc
        ## parallel multi-processing
        total_jobs = fa_count
        proc_list = []
        muscle_log_list = []

        ## determine how many jobs will be done at the same time
        if total_jobs <= number_proc:
            number_jobs = total_jobs
        else:
            number_jobs = number_proc

        ## spawn "number_jobs" jobs
        for job_index in xrange(number_jobs):

            ## call Popen to run MUSCLE from shell
            muscle_cmd,muscle_log_file = get_muscle_cmd_and_log(fa_files[job_index],outputdir)
            process = Popen(muscle_cmd,shell=True, stdout=PIPE, stderr=STDOUT)
            
            ## append process to proc_list; muscle_log to muscle_log_list
            proc_list.append(process)
            muscle_log_list.append(muscle_log_file)

        ## next job index
        job_index = number_jobs

        ## check if any job finished well, if so, spawn a new job to the finished process
        while (len(proc_list) > 0):

            ## loop proc_list
            for process_enum, process in enumerate(proc_list):                                                  
                                                                                                           
                ## check process status                                                                             
                process_status = process.poll()                                                                    
                                                                                                           
                ## if process is finished                                                                           
                if process_status is not None:                                                                     
                                                                                                           
                    # write muscle log file                                                                           
                    muscle_log_file = open(muscle_log_list[process_enum], 'wb')                                       
                    cmd_output = process.communicate()[0].strip()                                                  
                    muscle_log_file.write(cmd_output)                                                                  
                    muscle_log_file.close()                                                                            
                                                                                                           
                    # delete the process from the list                                                             
                    proc_list.pop(process_enum)                                                                 
                    muscle_log_list.pop(process_enum)                                                             
                                                                                                                   
                    # spawn another process until total jobs                                                       
                    if job_index < total_jobs:                                                                     
                                                                                                                   
                        # call Popen to run external script
                        muscle_cmd, muscle_log_file = get_muscle_cmd_and_log(fa_files[job_index],outputdir)
                        process = Popen(muscle_cmd,shell=True,stdout=PIPE, stderr = STDOUT)

                        # append process and samlog_path to list
                        proc_list.append(process)
                        muscle_log_list.append(muscle_log_file)

                        # job index + 1
                        job_index += 1

            # sleep a second
            time.sleep(1)

## ===========
## Main function
## ===========
def main(argv=None):
    # try to get arguments
    try:
        if argv is None:
            #argv = sys.argv[1:]
            args = parser.parse_args()
            print "arguments are {}".format(args)

            if args.quiet:
                print "running quietly"
            else:
                print "running verbosely"
                print "output directory is: {}".format(args.outputdir)
                if args.fastafiles:
                        if len(args.fastafiles) is 1:
                            print "input fasta file is: {}".format(args.fastafiles)
                        else:
                            print "input fasta files are: {}".format(args.fastafiles)
                if args.inputdirs:
                        if len(args.inputdirs) is 1:
                            print "input directory is: {}".format(args.inputdirs)
                        else:
                            print "input directories are: {}".format(args.inputdirs)
                if args.num_proc:
                    print "numebr of processes is: {}".format(args.num_proc)
            ## try to parse the arguments
            #option_map = muscle_parse_option(argv,version,help_message)
            start_time = datetime.now()
            sys.stderr.write("\n===============================================================================\n")
            sys.stderr.write("muscle starts running\n")
            muscleAlign(args)
            finish_time = datetime.now()
            duration = finish_time - start_time
            sys.stderr.write("Total elapsed time is %s [seconds]\n" %duration)
            sys.stderr.write("\n===============================================================================\n")
    

    ## Error handling
    except jj_utils.Usage, err:
        sys.stderr.write("%s: %s \n" %(os.path.basename(sys.argv[0]),str(err.msg)))
        sys.stderr.write("for help use -h/--help")
        return 2


if __name__ == "__main__":
    sys.exit(main())
