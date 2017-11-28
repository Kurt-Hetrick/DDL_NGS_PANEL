#$ -S /bin/bash
#$ -q cgc.q
#$ -cwd
#$ -V

umask 0007
umask
set

CORE_PATH="/isilon/ddl/SS_0500181/"
PICARD_DIR="/isilon/isilon-cgc/ddl_programs/picard-tools-1.103/"

PROJECT=$1
SAMPLE_NUMBER=$2
SM_TAG=$3

java -jar $PICARD_DIR/RevertSam.jar \
INPUT=$CORE_PATH/$PROJECT/$SM_TAG/BAM/$SM_TAG".bam" \
OUTPUT=/dev/stdout \
SORT_ORDER=queryname \
COMPRESSION_LEVEL=0 \
VALIDATION_STRINGENCY=SILENT | java -jar $PICARD_DIR/SamToFastq.jar \
INPUT=/dev/stdin \
FASTQ=$CORE_PATH/FASTQ_DEPO/$SM_TAG"_"$SAMPLE_NUMBER"_L001_R1_001.fastq.gz" \
SECOND_END_FASTQ=$CORE_PATH/FASTQ_DEPO/$SM_TAG"_"$SAMPLE_NUMBER"_L001_R2_001.fastq.gz" \
VALIDATION_STRINGENCY=SILENT
