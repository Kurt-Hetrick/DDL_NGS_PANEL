#!/bin/bash

CORE_PATH="/mnt/clinical/ddl/NGS/Panel_Data/"

BATCH=$1

mkdir -p $CORE_PATH/$BATCH/LOGS

sed 's/\r//g' $CORE_PATH/$BATCH/FASTQ/$BATCH".csv" \
| awk 'BEGIN {FS=","} NR>1 \
{print "qsub","-N",$1"_"$8"_CT12.001a.DDL.ClusterPipeline.v171201.sh","-o","/mnt/clinical/ddl/NGS/Panel_Data/"$1"/LOGS/"$1"_"$8".log",\
"-e","/mnt/clinical/ddl/NGS/Panel_Data/"$1"/LOGS/"$1"_"$8".log",\
"CT12.001a.DDL.ClusterPipeline.v171201.sh",\
$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$18,$19,$20"\n""sleep 3s"}'
