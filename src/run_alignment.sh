#!/bin/bash

#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=1
#PBS -q large
#PBS -N master_bbmap
#PBS -e /chongle/shared/work/pacbio_test/02_test_mock_community/01_PacBio_correction/master_bbmap.err
#PBS -o /chongle/shared/work/pacbio_test/02_test_mock_community/01_PacBio_correction/master_bbmap.out

cd /chongle/shared/work/pacbio_test/02_test_mock_community/01_PacBio_correction
echo Starting Time is $(date)

## total number of PacBio sequences
total_seqs=`grep '^>' /chongle/shared/database/03_PacBio/MockCommunity/PacBio/filtered_subreads/filtered_subreads.fasta -c`

## starting index of the PacBio sequences
seq_start=180000 # 0-based index of the PacBio sequence to align to

## total number of chunks to split the alignments to
num_chunks=$((($total_seqs - $seq_start) / 2000 ))

echo -e "\n"
echo total number of sequences is $total_seqs 
echo skip $seq_start sequences 
echo number of chunk jobs to submit: $num_chunks
echo -e "\n"

#underline=$'\137'

for i in `seq 0 $num_chunks`;
do
	running_jobs=`qstat -u cjg | grep R -c`
	queueing_jobs=`qstat -u cjg | grep R -c`
	total_jobs=$(($running_jobs + $queueing_jobs))
	until [ $total_jobs -lt 500 ]; do
		sleep 10m
		running_jobs=`qstat -u cjg | grep R -c`
		queueing_jobs=`qstat -u cjg | grep R -c`
		total_jobs=$(($running_jobs + $queueing_jobs))
	done
	echo chunk $i running
	skip_seqs=$(($seq_start + $i * 2000))
	start_seq=$(( $skip_seqs + 1))
	end_seq=$(( $skip_seqs + 2000 ))
	if [ $(($i % 2)) == 0 ]; then
		queue=medium
	else
		queue=large
	fi
	echo "python /chongle/shared/software/metalrec/src/check_PBReads.py -f /chongle/shared/database/03_PacBio/MockCommunity/PacBio/filtered_subreads/filtered_subreads.fasta-se /chongle/shared/work/pacbio_test/02_test_mock_community/illumina_subsample/flash_merge_sample_0.05.fastq -d /chongle/shared/work/pacbio_test/02_test_mock_community/01_PacBio_correction/mapping_${start_seq}_${end_seq} -k ${skip_seqs} -e $end_seq -q $queue > /chongle/shared/work/pacbio_test/02_test_mock_community/01_PacBio_correction/mapping_${start_seq}_${end_seq}.log"
	python /chongle/shared/software/metalrec/src/check_PBReads.py -f /chongle/shared/database/03_PacBio/MockCommunity/PacBio/filtered_subreads/filtered_subreads.fasta -se /chongle/shared/work/pacbio_test/02_test_mock_community/illumina_subsample/flash_merge_sample_0.05.fastq -d /chongle/shared/work/pacbio_test/02_test_mock_community/01_PacBio_correction/mapping_${start_seq}_${end_seq} -k ${skip_seqs} -e $end_seq -q $queue > /chongle/shared/work/pacbio_test/02_test_mock_community/01_PacBio_correction/mapping_${start_seq}_${end_seq}.log
done

echo Ending Time is $(date)
