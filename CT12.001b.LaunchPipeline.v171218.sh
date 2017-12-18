#!/bin/bash

CORE_PATH="/mnt/clinical/ddl/NGS/Panel_Data/"

BATCH=$1
QUEUE=$2

# If there are two arguments use the 2nd argument to define QUEUE
# If there is only one argument then cgc.q is the QUEUE

if [[ -z "$QUEUE" ]]
	then
	QUEUE=cgc.q
else
	QUEUE=$2
fi

# Make a directory for sge stdout/stderr logs

mkdir -p $CORE_PATH/$BATCH/LOGS

# make an array for each sample with information needed for pipeline input obtained from the sample sheet
# add a end of file is not present
# remove carriage returns if not present
# remove blank lines if present
# remove lines that only have whitespace

CREATE_SAMPLE_INFO_ARRAY ()
{
	SAMPLE_INFO_ARRAY=(`awk 1 $CORE_PATH/$BATCH/FASTQ/$BATCH".csv" | sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' | awk 'BEGIN{FS=","} $8=="'$SAMPLE'" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$18,$19,$20}' | sort | uniq`)
}

# I'd like to get away from explicitly calling the pipeline script here, but if they insist on the version being in the name, I can't think of any way around it.
# setting the priority to -1...gives 1 spot to move up if needed. but nobody else should be submitting at -1 afaik.
# previously priority was undefined so it was 0

PRINT_QSUB ()
{
	echo \
	qsub \
	-q $QUEUE \
	-p -1 \
	-N ${SAMPLE_INFO_ARRAY[0]}_${SAMPLE_INFO_ARRAY[7]}_CT12.001a.DDL.ClusterPipeline.v171201.sh \
	-j y \
	-o $CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/LOGS/${SAMPLE_INFO_ARRAY[0]}_${SAMPLE_INFO_ARRAY[7]}.log \
	CT12.001a.DDL.ClusterPipeline.v171201.sh \
	${SAMPLE_INFO_ARRAY[0]} \
	${SAMPLE_INFO_ARRAY[1]} \
	${SAMPLE_INFO_ARRAY[2]} \
	${SAMPLE_INFO_ARRAY[3]} \
	${SAMPLE_INFO_ARRAY[4]} \
	${SAMPLE_INFO_ARRAY[5]} \
	${SAMPLE_INFO_ARRAY[6]} \
	${SAMPLE_INFO_ARRAY[7]} \
	${SAMPLE_INFO_ARRAY[8]} \
	${SAMPLE_INFO_ARRAY[9]} \
	${SAMPLE_INFO_ARRAY[10]} \
	${SAMPLE_INFO_ARRAY[11]} \
	${SAMPLE_INFO_ARRAY[12]}
}

for SAMPLE in $(awk 1 $CORE_PATH/$BATCH/FASTQ/$BATCH".csv" | sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' | awk 'BEGIN {FS=","} NR>1 {print $8}' | sort | uniq );
do
	CREATE_SAMPLE_INFO_ARRAY
	PRINT_QSUB
	echo sleep 0.1s
done
