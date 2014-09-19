#!/bin/bash

#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=1
#PBS -q large
#PBS -N master_correct
#PBS -e /chongle/shared/work/pacbio_test/02_test_mock_community/01_PacBio_correction/master_correct.err
#PBS -o /chongle/shared/work/pacbio_test/02_test_mock_community/01_PacBio_correction/master_correct.out

cd /chongle/shared/work/pacbio_test/02_test_mock_community/01_PacBio_correction
echo Starting Time is $(date)

## total number of PacBio sequences
total_dirs=`find /chongle/shared/work/pacbio_test/02_test_mock_community/01_PacBio_correction/ -maxdepth 1 -mindepth 1 -type d`

echo -e "\n"

#underline=$'\137'

for dir in $total_dirs;
do
	total_jobs=`qstat -u cjg | grep 'R\|Q' -c`
	until [ $total_jobs -lt 2000 ]; do
		echo sleep 2m ...
		sleep 2m
		total_jobs=`qstat -u cjg | grep 'R\|Q' -c`
	done
	mediumQ=`qstat -u cjg | grep medium | grep Q -c`
	largeQ=`qstat -u cjg | grep large | grep Q -c`
	if [ $(( $largeQ - 2* $mediumQ )) -gt  0 ]; then
		queue=medium
	else
		queue=large
	fi
	#queue=medium
	echo $dir running on $queue queue
	echo "python /chongle/shared/software/metalrec/src/correct_PBReads.py -i $dir -q $queue -t 01:30:00"
	python /chongle/shared/software/metalrec/src/correct_PBReads.py -i $dir -q $queue -t 01:30:00
done

echo Ending Time is $(date)
