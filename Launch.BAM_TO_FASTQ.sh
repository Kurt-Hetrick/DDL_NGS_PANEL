#!/bin/bash

CORE_PATH="/isilon/ddl/SS_0500181/"

BATCH=$1

mkdir -p $CORE_PATH/$BATCH/LOGS

sed 's/\r//g' $CORE_PATH/$BATCH/FASTQ/$BATCH".csv" \
| awk 'BEGIN {FS=","} NR>1 \
{print "qsub","-N",$1"_"$8"_BAM_TO_FASTQ.sh","-o","/isilon/ddl/SS_0500181/"$1"/LOGS/"$1"_"$8".BAM_TO_FASTQ.log",\
"-e","/isilon/ddl/SS_0500181/"$1"/LOGS/"$1"_"$8".BAM_TO_FASTQ.log",\
"BAM_TO_FASTQ.sh",\
$1,$3,$8"\n""sleep 3s"}'
