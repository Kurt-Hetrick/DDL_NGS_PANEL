#$ -S /bin/bash
#$ -q cgc.q
#$ -cwd
#$ -V

umask 0007
umask
set

PROJECT=$1
FLOWCELL=$2
SAMPLE_NUMBER=$3
INDEX=$4
PLATFORM=$5
LIBRARY_NAME=$6
RUN_DATE=$7
SM_TAG=$8
CENTER=$9
DESCRIPTION=${10}
BAIT_BED=${11}
ROI_BED=${12}
ALL_TARGET_BED=${13}

PLATFORM_UNIT=$FLOWCELL"_"$SAMPLE_NUMBER"_"$INDEX

CORE_PATH="/isilon/ddl/SS_0500181/"
PIPE_FILES="/usr/local/sandbox/ddl_pipeline_files/"
BWA_DIR="/isilon/isilon-cgc/ddl_programs/bwa-0.7.5a/"
PICARD_DIR="/isilon/isilon-cgc/ddl_programs/picard-tools-1.103/"
GATK_DIR="/isilon/isilon-cgc/ddl_programs/GenomeAnalysisTK-2.7-4-g6f46d11/"
VERIFY_DIR="/isilon/isilon-cgc/ddl_programs/verifyBamID_20120620/bin/"
TABIX_DIR="/isilon/isilon-cgc/ddl_programs/tabix-0.2.6/"
BEDTOOLS_DIR="/isilon/isilon-cgc/ddl_programs/BEDTools-Version-2.16.2/bin/"
CIDRSEQSUITE_DIR="/isilon/isilon-cgc/ddl_programs/"
SAMTOOLS_DIR="/isilon/isilon-cgc/ddl_programs/samtools-0.1.18/"
REF_GENOME="human_g1k_v37_decoy.fasta"
BED_DIR=$CORE_PATH"/BED/"
GENE_LIST="RefSeqGene.GRCh37.Ready.txt"
BARCODE_BED="BarcodeSNPs.NGS1.v1.bed"
KNOWN_INDEL_1="Mills_and_1000G_gold_standard.indels.b37.vcf"
KNOWN_INDEL_2="1000G_phase1.indels.b37.vcf"
DBSNP="dbsnp_138.b37.vcf"
VERIFY_VCF="Omni25_genotypes_1525_samples_v2.b37.PASS.ALL.sites.vcf"

mkdir -p $CORE_PATH/$PROJECT/$SM_TAG/BAM
mkdir -p $CORE_PATH/$PROJECT/$SM_TAG/HC_BAM
mkdir -p $CORE_PATH/$PROJECT/$SM_TAG/ANALYSIS
mkdir -p $CORE_PATH/$PROJECT/$SM_TAG/Reports/{Alignment_Summary,ANNOVAR,Picard_Duplicates,QC_Reports,TiTv,verifyBamID}
mkdir -p $CORE_PATH/$PROJECT/$SM_TAG/Reports/BaseCall_Qscore_Distribution/{Metrics,PDF}
mkdir -p $CORE_PATH/$PROJECT/$SM_TAG/Reports/Count_Covariates/{CSV,PDF}
mkdir -p $CORE_PATH/$PROJECT/$SM_TAG/Reports/GC_Bias/{Metrics,PDF,Summary}
mkdir -p $CORE_PATH/$PROJECT/$SM_TAG/Reports/Genes_Coverage/{ROI,ALL_TARGET}
mkdir -p $CORE_PATH/$PROJECT/$SM_TAG/Reports/HybSelection/Per_Target_Coverage
mkdir -p $CORE_PATH/$PROJECT/$SM_TAG/Reports/Insert_Size/{Metrics,PDF}
mkdir -p $CORE_PATH/$PROJECT/$SM_TAG/Reports/Local_Realignment_Intervals
mkdir -p $CORE_PATH/$PROJECT/$SM_TAG/Reports/Mean_Quality_ByCycle/{Metrics,PDF}
mkdir -p $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/{RAW_BAIT,SNV_FILTERED_BAIT,SNV_FILTERED_ROI,INDEL_FILTERED_BAIT,INDEL_FILTERED_ROI,BARCODE}
mkdir -p $CORE_PATH/$PROJECT/$SM_TAG/temp

# -----Alignment and BAM post-processing-----

# --Bwa mem, trying to get all the reads from the fastq file out--
# --pipe to MergeSamFiles to sort and write a bam file.--

$BWA_DIR/bwa mem \
-M \
-T 0 \
$PIPE_FILES/$REF_GENOME \
$CORE_PATH/$PROJECT/FASTQ/$SM_TAG"_"$SAMPLE_NUMBER"_L001_R1_001.fastq.gz" \
$CORE_PATH/$PROJECT/FASTQ/$SM_TAG"_"$SAMPLE_NUMBER"_L001_R2_001.fastq.gz" \
-R "@RG\tID:$FLOWCELL\tPU:$PLATFORM_UNIT\tPL:$PLATFORM\tLB:$LIBRARY_NAME\tDT:$RUN_DATE\tSM:$SM_TAG\tCN:$CENTER\tDS:$DESCRIPTION" \
| java -jar $PICARD_DIR/MergeSamFiles.jar \
INPUT=/dev/stdin \
OUTPUT=$CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".original.bam" \
VALIDATION_STRINGENCY=SILENT \
SORT_ORDER=coordinate \
USE_THREADING=true \
CREATE_INDEX=true

## --Mark Duplicates with Picard, write a duplicate report

java -jar $PICARD_DIR/MarkDuplicates.jar \
INPUT=$CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".original.bam" \
OUTPUT=$CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".dup.bam" \
VALIDATION_STRINGENCY=SILENT \
METRICS_FILE=$CORE_PATH/$PROJECT/$SM_TAG/Reports/Picard_Duplicates/$SM_TAG".picard.duplicates.txt" \
CREATE_INDEX=true

## --Realigner Target Creator only on the baits, turn off downsampling

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-I $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".dup.bam" \
-R $PIPE_FILES/$REF_GENOME \
-L $BED_DIR/$BAIT_BED \
-known $PIPE_FILES/$KNOWN_INDEL_1 \
-known $PIPE_FILES/$KNOWN_INDEL_2 \
-dt NONE \
-o $CORE_PATH/$PROJECT/$SM_TAG/Reports/Local_Realignment_Intervals/$SM_TAG".intervals"

## --Local realignment turn off downsampling

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T IndelRealigner \
-I $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".dup.bam" \
-R $PIPE_FILES/$REF_GENOME \
-known $PIPE_FILES/$KNOWN_INDEL_1 \
-known $PIPE_FILES/$KNOWN_INDEL_2 \
-targetIntervals $CORE_PATH/$PROJECT/$SM_TAG/Reports/Local_Realignment_Intervals/$SM_TAG".intervals" \
-dt NONE \
-o $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".realign.bam"

## --BQSR using data only from the baited intervals, turn off downsampling--

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-I $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".realign.bam" \
-R $PIPE_FILES/$REF_GENOME \
-knownSites $PIPE_FILES/$KNOWN_INDEL_1 \
-knownSites $PIPE_FILES/$KNOWN_INDEL_2 \
-knownSites $PIPE_FILES/$DBSNP \
-dt NONE \
-L $BED_DIR/$BAIT_BED \
-o $CORE_PATH/$PROJECT/$SM_TAG/Reports/Count_Covariates/CSV/$SM_TAG".recal.table"

## --write out file with new scores, retain old scores, no downsampling

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T PrintReads \
-R $PIPE_FILES/$REF_GENOME \
-I $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".realign.bam" \
-BQSR $CORE_PATH/$PROJECT/$SM_TAG/Reports/Count_Covariates/CSV/$SM_TAG".recal.table" \
-dt NONE \
-EOQ \
-o $CORE_PATH/$PROJECT/$SM_TAG/BAM/$SM_TAG".bam"

# -----Reports Section-----

## --Alignment Summary Metrics--

java -jar $PICARD_DIR/CollectAlignmentSummaryMetrics.jar \
INPUT=$CORE_PATH/$PROJECT/$SM_TAG/BAM/$SM_TAG".bam" \
OUTPUT=$CORE_PATH/$PROJECT/$SM_TAG/Reports/Alignment_Summary/$SM_TAG"_alignment_summary_metrics.txt" \
R=$PIPE_FILES/$REF_GENOME \
VALIDATION_STRINGENCY=SILENT

## --Base Call Quality Score Distribution--

java -jar $PICARD_DIR/QualityScoreDistribution.jar \
INPUT=$CORE_PATH/$PROJECT/$SM_TAG/BAM/$SM_TAG".bam" \
OUTPUT=$CORE_PATH/$PROJECT/$SM_TAG/Reports/BaseCall_Qscore_Distribution/Metrics/$SM_TAG"_quality_score_distribution.txt" \
CHART=$CORE_PATH/$PROJECT/$SM_TAG/Reports/BaseCall_Qscore_Distribution/PDF/$SM_TAG"_quality_score_distribution.pdf" \
R=$PIPE_FILES/$REF_GENOME \
VALIDATION_STRINGENCY=SILENT \
INCLUDE_NO_CALLS=true

## --GC Bias Metrics--

java -jar $PICARD_DIR/CollectGcBiasMetrics.jar \
INPUT=$CORE_PATH/$PROJECT/$SM_TAG/BAM/$SM_TAG".bam" \
OUTPUT=$CORE_PATH/$PROJECT/$SM_TAG/Reports/GC_Bias/Metrics/$SM_TAG"_gc_bias_metrics.txt" \
CHART_OUTPUT=$CORE_PATH/$PROJECT/$SM_TAG/Reports/GC_Bias/PDF/$SM_TAG"_gc_bias_metrics.pdf" \
SUMMARY_OUTPUT=$CORE_PATH/$PROJECT/$SM_TAG/Reports/GC_Bias/Summary/$SM_TAG"_gc_bias_summary.txt" \
R=$PIPE_FILES/$REF_GENOME \
VALIDATION_STRINGENCY=SILENT

## --Insert Size--

java -jar $PICARD_DIR/CollectInsertSizeMetrics.jar \
INPUT=$CORE_PATH/$PROJECT/$SM_TAG/BAM/$SM_TAG".bam" \
OUTPUT=$CORE_PATH/$PROJECT/$SM_TAG/Reports/Insert_Size/Metrics/$SM_TAG"_insert_size_metrics.txt" \
H=$CORE_PATH/$PROJECT/$SM_TAG/Reports/Insert_Size/PDF/$SM_TAG"_insert_size_metrics_histogram.pdf" \
R=$PIPE_FILES/$REF_GENOME \
VALIDATION_STRINGENCY=SILENT

## --Mean Quality By Cycle--

java -jar $PICARD_DIR/MeanQualityByCycle.jar \
INPUT=$CORE_PATH/$PROJECT/$SM_TAG/BAM/$SM_TAG".bam" \
OUTPUT=$CORE_PATH/$PROJECT/$SM_TAG/Reports/Mean_Quality_ByCycle/Metrics/$SM_TAG"_mean_quality_by_cycle.txt" \
CHART=$CORE_PATH/$PROJECT/$SM_TAG/Reports/Mean_Quality_ByCycle/PDF/$SM_TAG"_mean_quality_by_cycle_chart.pdf" \
R=$PIPE_FILES/$REF_GENOME \
VALIDATION_STRINGENCY=SILENT

## --Depth of Coverage for the regions of interest--

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R $PIPE_FILES/$REF_GENOME \
-I $CORE_PATH/$PROJECT/$SM_TAG/BAM/$SM_TAG".bam" \
-geneList:REFSEQ $PIPE_FILES/$GENE_LIST \
-L $BED_DIR/$ROI_BED \
-mmq 17 \
-mbq 20 \
--outputFormat table \
-o $CORE_PATH/$PROJECT/$SM_TAG/Reports/Genes_Coverage/ROI/$SM_TAG"_"$ROI_BED \
-ct 30 \
-ct 50

## Removing carriage returns in the ROI bed file.

sed 's/\r//g' $BED_DIR/$ROI_BED > $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG"_"$ROI_BED

## Okay, reformatting per-base interval FOR ROI

awk 'BEGIN {OFS="\t"} NR>1 {split($1,BASE,":"); print BASE[1],BASE[2]-1,BASE[2],$2}' \
$CORE_PATH/$PROJECT/$SM_TAG/Reports/Genes_Coverage/ROI/$SM_TAG"_"$ROI_BED \
| $BEDTOOLS_DIR/bedtools intersect -a $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG"_"$ROI_BED -b /dev/stdin -wb \
| awk '{OFS="\t"} BEGIN {print "CHROM""\t""POS""\t""ROI""\t""STRAND""\t""DEPTH"} {print $1,$3,$4,$6,$10}' \
>| $CORE_PATH/$PROJECT/$SM_TAG/Reports/Genes_Coverage/ROI/$SM_TAG"_"$ROI_BED".PER.BASE.REPORT.txt"

cp -rvf $CORE_PATH/$PROJECT/$SM_TAG/Reports/Genes_Coverage/ROI/$SM_TAG"_"$ROI_BED".PER.BASE.REPORT.txt" \
$CORE_PATH/$PROJECT/$SM_TAG/ANALYSIS/$SM_TAG"_"$ROI_BED".PER.BASE.REPORT.txt"

## Now reformatting sample interval summary FOR ROI

awk 'NR>1' $CORE_PATH/$PROJECT/$SM_TAG/Reports/Genes_Coverage/ROI/$SM_TAG"_"$ROI_BED".sample_interval_summary" \
| sed 's/:/\t/g' \
| awk 'BEGIN {OFS="\t"} $2!~"-" {print $1,$2"-"$2,$3,$4,$5,$6,$7,$8,$9,$10,$11} $2~"-" {print $0}' \
| awk 'BEGIN {OFS="\t"} {gsub (/-/,"\t"); print $1,$2,$3,$4,$5,$8,$9,$10,$11,$12}' \
| $BEDTOOLS_DIR/bedtools intersect -a $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG"_"$ROI_BED -b /dev/stdin -wb \
| awk '{OFS="\t"} BEGIN {print "CHROM""\t""START""\t""END""\t""ROI""\t""STRAND""\t""TOTAL_CVG""\t""AVG_CVG""\t""Q1_CVG""\t""MEDIAN_CVG""\t""Q3_CVG""\t""PCT_30x""\t""PCT_50x"} \
{print $7,$8,$9,$4,$6,$10,$11,$12,$13,$14,$15,$16}' \
>| $CORE_PATH/$PROJECT/$SM_TAG/Reports/Genes_Coverage/ROI/$SM_TAG"_"$ROI_BED".PER.INTERVAL.SUMMARY.txt"

cp -rvf $CORE_PATH/$PROJECT/$SM_TAG/Reports/Genes_Coverage/ROI/$SM_TAG"_"$ROI_BED".PER.INTERVAL.SUMMARY.txt" \
$CORE_PATH/$PROJECT/$SM_TAG/ANALYSIS/$SM_TAG"_"$ROI_BED".PER.INTERVAL.SUMMARY.txt"

## Same as above but only outputting those records where the PCT_50x is less than 100 FOR ROI

awk 'NR>1' $CORE_PATH/$PROJECT/$SM_TAG/Reports/Genes_Coverage/ROI/$SM_TAG"_"$ROI_BED".sample_interval_summary" \
| sed 's/:/\t/g' \
| awk 'BEGIN {OFS="\t"} $2!~"-" {print $1,$2"-"$2,$3,$4,$5,$6,$7,$8,$9,$10,$11} $2~"-" {print $0}' \
| awk 'BEGIN {OFS="\t"} {gsub (/-/,"\t"); print $1,$2,$3,$4,$5,$8,$9,$10,$11,$12}' \
| $BEDTOOLS_DIR/bedtools intersect -a $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG"_"$ROI_BED -b /dev/stdin -wb \
| awk '{OFS="\t"} BEGIN {print "CHROM""\t""START""\t""END""\t""ROI""\t""STRAND""\t""TOTAL_CVG""\t""AVG_CVG""\t""Q1_CVG""\t""MEDIAN_CVG""\t""Q3_CVG""\t""PCT_30x""\t""PCT_50x"} \
$16<100 \
{print $7,$8,$9,$4,$6,$10,$11,$12,$13,$14,$15,$16}' \
>| $CORE_PATH/$PROJECT/$SM_TAG/Reports/Genes_Coverage/ROI/$SM_TAG"_"$ROI_BED".PER.INTERVAL.SUMMARY.LT.50x.txt"

cp -rvf $CORE_PATH/$PROJECT/$SM_TAG/Reports/Genes_Coverage/ROI/$SM_TAG"_"$ROI_BED".PER.INTERVAL.SUMMARY.LT.50x.txt" \
$CORE_PATH/$PROJECT/$SM_TAG/ANALYSIS/$SM_TAG"_"$ROI_BED".PER.INTERVAL.SUMMARY.LT.50x.txt"

## --Depth of Coverage for the all targets bed file--

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R $PIPE_FILES/$REF_GENOME \
-I $CORE_PATH/$PROJECT/$SM_TAG/BAM/$SM_TAG".bam" \
-geneList:REFSEQ $PIPE_FILES/$GENE_LIST \
-L $BED_DIR/$ALL_TARGET_BED \
-mmq 17 \
-mbq 20 \
--outputFormat table \
-o $CORE_PATH/$PROJECT/$SM_TAG/Reports/Genes_Coverage/ALL_TARGET/$SM_TAG"_"$ALL_TARGET_BED \
-ct 30 \
-ct 50

# Removing carriage returns in the all target bed file

sed 's/\r//g' $BED_DIR/$ALL_TARGET_BED > $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG"_"$ALL_TARGET_BED

## Okay, reformatting per-base interval FOR ALL TARGETS

awk 'BEGIN {OFS="\t"} NR>1 {split($1,BASE,":"); print BASE[1],BASE[2]-1,BASE[2],$2}' \
$CORE_PATH/$PROJECT/$SM_TAG/Reports/Genes_Coverage/ALL_TARGET/$SM_TAG"_"$ALL_TARGET_BED \
| $BEDTOOLS_DIR/bedtools intersect -a $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG"_"$ALL_TARGET_BED -b /dev/stdin -wb \
| awk '{OFS="\t"} BEGIN {print "CHROM""\t""POS""\t""ROI""\t""STRAND""\t""DEPTH"} {print $1,$3,$4,$6,$10}' \
>| $CORE_PATH/$PROJECT/$SM_TAG/Reports/Genes_Coverage/ALL_TARGET/$SM_TAG"_"$ALL_TARGET_BED".PER.BASE.REPORT.txt"

cp -rvf $CORE_PATH/$PROJECT/$SM_TAG/Reports/Genes_Coverage/ALL_TARGET/$SM_TAG"_"$ALL_TARGET_BED".PER.BASE.REPORT.txt" \
$CORE_PATH/$PROJECT/$SM_TAG/ANALYSIS/$SM_TAG"_"$ALL_TARGET_BED".PER.BASE.REPORT.txt"

## Now reformatting sample interval summary FOR ALL TARGETS

awk 'NR>1' $CORE_PATH/$PROJECT/$SM_TAG/Reports/Genes_Coverage/ALL_TARGET/$SM_TAG"_"$ALL_TARGET_BED".sample_interval_summary" \
| sed 's/:/\t/g' \
| awk 'BEGIN {OFS="\t"} $2!~"-" {print $1,$2"-"$2,$3,$4,$5,$6,$7,$8,$9,$10,$11} $2~"-" {print $0}' \
| awk 'BEGIN {OFS="\t"} {gsub (/-/,"\t"); print $1,$2,$3,$4,$5,$8,$9,$10,$11,$12}' \
| $BEDTOOLS_DIR/bedtools intersect -a $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG"_"$ALL_TARGET_BED -b /dev/stdin -wb \
| awk '{OFS="\t"} BEGIN {print "CHROM""\t""START""\t""END""\t""ROI""\t""STRAND""\t""TOTAL_CVG""\t""AVG_CVG""\t""Q1_CVG""\t""MEDIAN_CVG""\t""Q3_CVG""\t""PCT_30x""\t""PCT_50x"} \
{print $7,$8,$9,$4,$6,$10,$11,$12,$13,$14,$15,$16}' \
>| $CORE_PATH/$PROJECT/$SM_TAG/Reports/Genes_Coverage/ALL_TARGET/$SM_TAG"_"$ALL_TARGET_BED".PER.INTERVAL.SUMMARY.txt"

cp -rvf $CORE_PATH/$PROJECT/$SM_TAG/Reports/Genes_Coverage/ALL_TARGET/$SM_TAG"_"$ALL_TARGET_BED".PER.INTERVAL.SUMMARY.txt" \
$CORE_PATH/$PROJECT/$SM_TAG/ANALYSIS/$SM_TAG"_"$ALL_TARGET_BED".PER.INTERVAL.SUMMARY.txt"

## Same as above but only outputting those records where the PCT_50x is less than 100 FOR ALL TARGETS

awk 'NR>1' $CORE_PATH/$PROJECT/$SM_TAG/Reports/Genes_Coverage/ALL_TARGET/$SM_TAG"_"$ALL_TARGET_BED".sample_interval_summary" \
| sed 's/:/\t/g' \
| awk 'BEGIN {OFS="\t"} $2!~"-" {print $1,$2"-"$2,$3,$4,$5,$6,$7,$8,$9,$10,$11} $2~"-" {print $0}' \
| awk 'BEGIN {OFS="\t"} {gsub (/-/,"\t"); print $1,$2,$3,$4,$5,$8,$9,$10,$11,$12}' \
| $BEDTOOLS_DIR/bedtools intersect -a $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG"_"$ALL_TARGET_BED -b /dev/stdin -wb \
| awk '{OFS="\t"} BEGIN {print "CHROM""\t""START""\t""END""\t""ROI""\t""STRAND""\t""TOTAL_CVG""\t""AVG_CVG""\t""Q1_CVG""\t""MEDIAN_CVG""\t""Q3_CVG""\t""PCT_30x""\t""PCT_50x"} \
$16<100 \
{print $7,$8,$9,$4,$6,$10,$11,$12,$13,$14,$15,$16}' \
>| $CORE_PATH/$PROJECT/$SM_TAG/Reports/Genes_Coverage/ALL_TARGET/$SM_TAG"_"$ALL_TARGET_BED".PER.INTERVAL.SUMMARY.LT.50x.txt"

cp -rvf $CORE_PATH/$PROJECT/$SM_TAG/Reports/Genes_Coverage/ALL_TARGET/$SM_TAG"_"$ALL_TARGET_BED".PER.INTERVAL.SUMMARY.LT.50x.txt" \
$CORE_PATH/$PROJECT/$SM_TAG/ANALYSIS/$SM_TAG"_"$ALL_TARGET_BED".PER.INTERVAL.SUMMARY.LT.50x.txt"


## --Creating "all targets" and "on target" picard bed files for use in CalculateHS metrics--

( $SAMTOOLS_DIR/samtools view -H $CORE_PATH/$PROJECT/$SM_TAG/BAM/$SM_TAG".bam" | grep "@SQ" ; sed 's/\r//g' $BED_DIR/$ALL_TARGET_BED | awk 'NR>1 {print $1,($2+1),$3,"+",$1"_"($2+1)"_"$3}' | sed 's/ /\t/g') \
>| $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG"_"$ALL_TARGET_BED".picard.bed"

( $SAMTOOLS_DIR/samtools view -H $CORE_PATH/$PROJECT/$SM_TAG/BAM/$SM_TAG".bam" | grep "@SQ" ; sed 's/\r//g' $BED_DIR/$ROI_BED | awk 'NR>1 {print $1,($2+1),$3,"+",$1"_"($2+1)"_"$3}' | sed 's/ /\t/g') \
>| $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG"_"$ROI_BED".picard.bed"

## --Calculate HS metrics and also calculate a per target coverage report--

java -jar $PICARD_DIR/CalculateHsMetrics.jar \
INPUT=$CORE_PATH/$PROJECT/$SM_TAG/BAM/$SM_TAG".bam" \
OUTPUT=$CORE_PATH/$PROJECT/$SM_TAG/Reports/HybSelection/$SM_TAG"_hybridization_selection_metrics.txt" \
PER_TARGET_COVERAGE=$CORE_PATH/$PROJECT/$SM_TAG/Reports/HybSelection/Per_Target_Coverage/$SM_TAG"_per_target_coverage.txt" \
R=$PIPE_FILES/$REF_GENOME \
BI=$CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG"_"$ALL_TARGET_BED".picard.bed" \
TI=$CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG"_"$ROI_BED".picard.bed" \
VALIDATION_STRINGENCY=SILENT

## --Creating a VCF file to be used as the reference for verifyBamID--

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $PIPE_FILES/$REF_GENOME \
--variant $PIPE_FILES/$VERIFY_VCF \
-L $BED_DIR/$BAIT_BED \
-o $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".VerifyBamID.vcf"

## --Running verifyBamID--

$VERIFY_DIR/verifyBamID \
--bam $CORE_PATH/$PROJECT/$SM_TAG/BAM/$SM_TAG".bam" \
--vcf $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".VerifyBamID.vcf" \
--out $CORE_PATH/$PROJECT/$SM_TAG/Reports/verifyBamID/$SM_TAG \
--ignoreRG \
--precise \
--verbose \
--maxDepth 2000

## Eventually Below will get wrapped in the QC report

awk 'BEGIN {OFS="\t"} NR==1 {print $1,$4,$5,$6,$7,$8,$9,"DIFF_LK0_LK1"} \
NR==2 {print $1,$4,$5,$6,$7,$8,$9,$9-$8}' \
$CORE_PATH/$PROJECT/$SM_TAG/Reports/verifyBamID/$SM_TAG".selfSM" \
>| $CORE_PATH/$PROJECT/$SM_TAG/ANALYSIS/$SM_TAG".VerifyBamID.txt"

## --Generate post BQSR table--

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-I $CORE_PATH/$PROJECT/$SM_TAG/BAM/$SM_TAG".bam" \
-R $PIPE_FILES/$REF_GENOME \
-knownSites $PIPE_FILES/$KNOWN_INDEL_1 \
-knownSites $PIPE_FILES/$KNOWN_INDEL_2 \
-knownSites $PIPE_FILES/$DBSNP \
-dt NONE \
-BQSR $CORE_PATH/$PROJECT/$SM_TAG/Reports/Count_Covariates/CSV/$SM_TAG".recal.table" \
-o $CORE_PATH/$PROJECT/$SM_TAG/Reports/Count_Covariates/CSV/$SM_TAG".recal.post.table" \
-nct 8

## --Generate BQSR plots--

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T AnalyzeCovariates \
-R $PIPE_FILES/$REF_GENOME \
-before $CORE_PATH/$PROJECT/$SM_TAG/Reports/Count_Covariates/CSV/$SM_TAG".recal.table" \
-after $CORE_PATH/$PROJECT/$SM_TAG/Reports/Count_Covariates/CSV/$SM_TAG".recal.post.table" \
-plots $CORE_PATH/$PROJECT/$SM_TAG/Reports/Count_Covariates/PDF/$SM_TAG".bqsr.pdf"


## -----Unified Genotyper Section-----

## --Unified Genotyper, Emit All Sites, SNV model only, Emit and make calls at QUAL>=0, min baseq at 20, no downsampling--
## --Call on padded bait--

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R $PIPE_FILES/$REF_GENOME \
--input_file $CORE_PATH/$PROJECT/$SM_TAG/BAM/$SM_TAG".bam" \
--dbsnp $PIPE_FILES/$DBSNP \
-L $BED_DIR/$BAIT_BED \
-glm SNP \
--output_mode EMIT_ALL_SITES \
-stand_emit_conf 0 \
-stand_call_conf 0 \
-mbq 20 \
-A QualByDepth \
-A HaplotypeScore \
-A MappingQualityRankSumTest \
-A ReadPosRankSumTest \
-A FisherStrand \
-A GCContent \
-A AlleleBalanceBySample \
-A AlleleBalance \
-A VariantType \
-dt NONE \
-o $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/RAW_BAIT/$SM_TAG".UG.EAS.raw.OnBait.SNV.vcf"

## --Filtering "SNVs"--
## --Logging_level ERROR is used here--

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R $PIPE_FILES/$REF_GENOME \
--variant $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/RAW_BAIT/$SM_TAG".UG.EAS.raw.OnBait.SNV.vcf" \
--filterExpression 'QD < 2.0' \
--filterName 'QDfilter' \
--filterExpression 'ABHet > 0.80' \
--filterName 'ABfilter80' \
--filterExpression 'ABHet < 0.20' \
--filterName 'ABfilter20' \
--filterExpression 'QUAL < 30.0' \
--filterName 'QUALfilter' \
--filterExpression 'FS > 20.0' \
--filterName 'FSfilter' \
--logging_level ERROR \
-o $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/SNV_FILTERED_BAIT/$SM_TAG".UG.EAS.FILTERED.OnBait.SNV.vcf"

## --bgzip and tabix index "UG SNV", use tabix to filter on ROI--

$TABIX_DIR/bgzip -c $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/SNV_FILTERED_BAIT/$SM_TAG".UG.EAS.FILTERED.OnBait.SNV.vcf" \
>| $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/SNV_FILTERED_BAIT/$SM_TAG".UG.EAS.FILTERED.OnBait.SNV.vcf.gz"

$TABIX_DIR/tabix -f -p vcf $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/SNV_FILTERED_BAIT/$SM_TAG".UG.EAS.FILTERED.OnBait.SNV.vcf.gz"

$TABIX_DIR/tabix -h $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/SNV_FILTERED_BAIT/$SM_TAG".UG.EAS.FILTERED.OnBait.SNV.vcf.gz" \
-B $BED_DIR/$ROI_BED \
>| $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/SNV_FILTERED_ROI/$SM_TAG".UG.EAS.FILTERED.SNV."$ROI_BED".vcf"

echo Start ANNOVAR UG `date`

java -jar $CIDRSEQSUITE_DIR/CIDRSeqSuite.jar -pipeline -annovar_annotation \
$CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/SNV_FILTERED_ROI/$SM_TAG".UG.EAS.FILTERED.SNV."$ROI_BED".vcf" \
$CORE_PATH/$PROJECT/$SM_TAG/temp/ 1

echo End ANNOVAR UG `date`

cp -rvf $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".UG.EAS.FILTERED.SNV."$ROI_BED"_ANNOVAR_DICTIONARY_2013_02_21.csv" \
$CORE_PATH/$PROJECT/$SM_TAG/Reports/ANNOVAR/$SM_TAG".UG.EAS.FILTERED.SNV."$ROI_BED"_ANNOVAR_DICTIONARY_2013_02_21.csv"

cp -rvf $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".UG.EAS.FILTERED.SNV."$ROI_BED"_ANNOVAR_REPORT.txt" \
$CORE_PATH/$PROJECT/$SM_TAG/Reports/ANNOVAR/$SM_TAG".UG.EAS.FILTERED.SNV."$ROI_BED"_ANNOVAR_REPORT.txt"

cp -rvf $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".UG.EAS.FILTERED.SNV."$ROI_BED"_ANNOVAR_REPORT.txt" \
$CORE_PATH/$PROJECT/$SM_TAG/ANALYSIS/$SM_TAG".UG.EAS.FILTERED.SNV."$ROI_BED"_ANNOVAR_REPORT.txt"

cp -rvf $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/SNV_FILTERED_ROI/$SM_TAG".UG.EAS.FILTERED.SNV."$ROI_BED".vcf" \
$CORE_PATH/$PROJECT/$SM_TAG/ANALYSIS/$SM_TAG".UG.EAS.FILTERED.SNV."$ROI_BED".vcf"

## Extracting the 91 Barcode SNPs

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $PIPE_FILES/$REF_GENOME \
--variant $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/SNV_FILTERED_BAIT/$SM_TAG".UG.EAS.FILTERED.OnBait.SNV.vcf" \
-selectType SNP \
-selectType NO_VARIATION \
-L $BED_DIR/$BARCODE_BED \
-o $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/BARCODE/$SM_TAG".UG.EAS.FILTERED.BARCODE.SNV.vcf"

cp -rvf $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/BARCODE/$SM_TAG".UG.EAS.FILTERED.BARCODE.SNV.vcf" \
$CORE_PATH/$PROJECT/$SM_TAG/ANALYSIS/$SM_TAG".UG.EAS.FILTERED.BARCODE.SNV.vcf"

## -----Haplotype Caller-----

## Call on Bait (padded)

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R $PIPE_FILES/$REF_GENOME \
--input_file $CORE_PATH/$PROJECT/$SM_TAG/BAM/$SM_TAG".bam" \
--dbsnp $PIPE_FILES/$DBSNP \
-L $BED_DIR/$BAIT_BED \
-stand_emit_conf 0 \
-stand_call_conf 0 \
-A VariantType \
-A DepthPerSampleHC \
-A ClippingRankSumTest \
-A MappingQualityRankSumTest \
-A ReadPosRankSumTest \
-A FisherStrand \
-A GCContent \
-A AlleleBalanceBySample \
-A AlleleBalance \
-A QualByDepth \
-dt NONE \
-o $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".HC.raw.OnBait.vcf"

## --HC does not output AlleleBalance--
## --VariantAnnotator will add in the ABHet--

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T VariantAnnotator \
-R $PIPE_FILES/$REF_GENOME \
--input_file $CORE_PATH/$PROJECT/$SM_TAG/BAM/$SM_TAG".bam" \
--variant $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".HC.raw.OnBait.vcf" \
-L $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".HC.raw.OnBait.vcf" \
-A AlleleBalance \
-A VariantType \
-dt NONE \
-o $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".HC.raw.OnBait.ANNOTATED.vcf"

# --Extracting INDELs--

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $PIPE_FILES/$REF_GENOME \
--variant $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".HC.raw.OnBait.ANNOTATED.vcf" \
-selectType INDEL \
-o $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/RAW_BAIT/$SM_TAG".HC.raw.OnBait.INDEL.vcf"

## --Filtering "INDELs"--

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R $PIPE_FILES/$REF_GENOME \
--variant $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/RAW_BAIT/$SM_TAG".HC.raw.OnBait.INDEL.vcf" \
--filterExpression 'QD < 2.0 && QUAL < 1000.0' \
--filterName 'QD_QUAL_filter_INDEL' \
--filterExpression 'ReadPosRankSum < -20.0' \
--filterName 'ReadPosRankSum_INDEL' \
--filterExpression 'FS > 200.0' \
--filterName 'FSfilter_INDEL' \
--filterExpression 'QUAL < 30.0' \
--filterName 'QUALfilter_INDEL' \
--logging_level ERROR \
-o $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/INDEL_FILTERED_BAIT/$SM_TAG".HC.FILTERED.OnBait.INDEL.vcf"

## --bgzip and tabix index "HC indel calls", use tabix to filter on ROI--

$TABIX_DIR/bgzip -c $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/INDEL_FILTERED_BAIT/$SM_TAG".HC.FILTERED.OnBait.INDEL.vcf" \
>| $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/INDEL_FILTERED_BAIT/$SM_TAG".HC.FILTERED.OnBait.INDEL.vcf.gz"

$TABIX_DIR/tabix -f -p vcf $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/INDEL_FILTERED_BAIT/$SM_TAG".HC.FILTERED.OnBait.INDEL.vcf.gz"

$TABIX_DIR/tabix -h $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/INDEL_FILTERED_BAIT/$SM_TAG".HC.FILTERED.OnBait.INDEL.vcf.gz" \
-B $BED_DIR/$ROI_BED \
>| $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/INDEL_FILTERED_ROI/$SM_TAG".HC.FILTERED.INDEL."$ROI_BED".vcf"

## Extracting MNP, SYMBOLIC, MIXED

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $PIPE_FILES/$REF_GENOME \
--variant $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".HC.raw.OnBait.ANNOTATED.vcf" \
-selectType MIXED \
-selectType MNP \
-selectType SYMBOLIC \
-o $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/RAW_BAIT/$SM_TAG".HC.raw.OnBait.OTHER.vcf"

## --Filtering "MNP, SYMBOLIC, MIXED"--

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R $PIPE_FILES/$REF_GENOME \
--variant $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/RAW_BAIT/$SM_TAG".HC.raw.OnBait.OTHER.vcf" \
--filterExpression 'QD < 2.0 && QUAL < 1000.0' \
--filterName 'QD_QUAL_filter_INDEL' \
--filterExpression 'ReadPosRankSum < -20.0' \
--filterName 'ReadPosRankSum_INDEL' \
--filterExpression 'FS > 200.0' \
--filterName 'FSfilter_INDEL' \
--filterExpression 'QUAL < 30.0' \
--filterName 'QUALfilter_INDEL' \
--logging_level ERROR \
-o $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/INDEL_FILTERED_BAIT/$SM_TAG".HC.FILTERED.OnBait.OTHER.vcf"

## --bgzip and tabix index "HC OTHER calls", use tabix to filter on ROI--

$TABIX_DIR/bgzip -c $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/INDEL_FILTERED_BAIT/$SM_TAG".HC.FILTERED.OnBait.OTHER.vcf" \
>| $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/INDEL_FILTERED_BAIT/$SM_TAG".HC.FILTERED.OnBait.OTHER.vcf.gz"

$TABIX_DIR/tabix -f -p vcf $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/INDEL_FILTERED_BAIT/$SM_TAG".HC.FILTERED.OnBait.OTHER.vcf.gz"

awk 'NR==1 {print "'$TABIX_DIR'""/tabix -h",\
"'$CORE_PATH'""/""'$PROJECT'""/""'$SM_TAG'""/VARIANT_CALLS/INDEL_FILTERED_BAIT/""'$SM_TAG'"".HC.FILTERED.OnBait.OTHER.vcf.gz",\
$1":"$2"-"$3}' \
$CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG"_"$ROI_BED \
| bash >| $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/INDEL_FILTERED_ROI/$SM_TAG".HC.FILTERED.OTHER."$ROI_BED".vcf"

awk 'NR>1 {print "'$TABIX_DIR'""/tabix",\
"'$CORE_PATH'""/""'$PROJECT'""/""'$SM_TAG'""/VARIANT_CALLS/INDEL_FILTERED_BAIT/""'$SM_TAG'"".HC.FILTERED.OnBait.OTHER.vcf.gz",\
$1":"$2"-"$3}' \
$CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG"_"$ROI_BED \
| bash >> $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/INDEL_FILTERED_ROI/$SM_TAG".HC.FILTERED.OTHER."$ROI_BED".vcf"

cp -rvf $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/INDEL_FILTERED_ROI/$SM_TAG".HC.FILTERED.OTHER."$ROI_BED".vcf" \
$CORE_PATH/$PROJECT/$SM_TAG/ANALYSIS/$SM_TAG".HC.FILTERED.OTHER."$ROI_BED".vcf"

echo Start ANNOVAR HC `date`

java -jar $CIDRSEQSUITE_DIR/CIDRSeqSuite.jar -pipeline -annovar_annotation \
$CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/INDEL_FILTERED_ROI/$SM_TAG".HC.FILTERED.INDEL."$ROI_BED".vcf" \
$CORE_PATH/$PROJECT/$SM_TAG/temp/ 1

echo End ANNOVAR HC `date`

cp -rvf $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".HC.FILTERED.INDEL."$ROI_BED"_ANNOVAR_DICTIONARY_2013_02_21.csv" \
$CORE_PATH/$PROJECT/$SM_TAG/Reports/ANNOVAR/$SM_TAG".HC.FILTERED.INDEL."$ROI_BED"_ANNOVAR_DICTIONARY_2013_02_21.csv"

cp -rvf $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".HC.FILTERED.INDEL."$ROI_BED"_ANNOVAR_REPORT.txt" \
$CORE_PATH/$PROJECT/$SM_TAG/Reports/ANNOVAR/$SM_TAG".HC.FILTERED.INDEL."$ROI_BED"_ANNOVAR_REPORT.txt"

cp -rvf $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".HC.FILTERED.INDEL."$ROI_BED"_ANNOVAR_REPORT.txt" \
$CORE_PATH/$PROJECT/$SM_TAG/ANALYSIS/$SM_TAG".HC.FILTERED.INDEL."$ROI_BED"_ANNOVAR_REPORT.txt"

cp -rvf $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/INDEL_FILTERED_ROI/$SM_TAG".HC.FILTERED.INDEL."$ROI_BED".vcf" \
$CORE_PATH/$PROJECT/$SM_TAG/ANALYSIS/$SM_TAG".HC.FILTERED.INDEL."$ROI_BED".vcf"

## Create a bed file for the HC ROI indel calls with a 150bp pad from the leftmost coordinate.

grep -v "^#" $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/INDEL_FILTERED_ROI/$SM_TAG".HC.FILTERED.INDEL."$ROI_BED".vcf" \
| awk 'BEGIN {OFS="\t"} {print $1,$2-150,$2+150}' > $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".HC.FILTERED.INDEL."$ROI_BED".vcf.bed"

## Use the above bed file to write a HC bam file for ROI indels

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R $PIPE_FILES/$REF_GENOME \
--input_file $CORE_PATH/$PROJECT/$SM_TAG/BAM/$SM_TAG".bam" \
--dbsnp $PIPE_FILES/$DBSNP \
-L $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".HC.FILTERED.INDEL."$ROI_BED".vcf.bed" \
--bamOutput $CORE_PATH/$PROJECT/$SM_TAG/HC_BAM/$SM_TAG".HC.FILTERED."$ROI_BED".INDEL.bam" \
-o $CORE_PATH/$PROJECT/$SM_TAG/HC_BAM/$SM_TAG".HC.BAM."$ROI_BED".vcf"

## Create a vcf file for passing UG SNVs on bait

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $PIPE_FILES/$REF_GENOME \
--variant $CORE_PATH/$PROJECT/$SM_TAG/VARIANT_CALLS/SNV_FILTERED_BAIT/$SM_TAG".UG.EAS.FILTERED.OnBait.SNV.vcf" \
-selectType SNP \
-ef \
-env \
-o $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".UG.EAS.FILTERED.OnBait.TiTv.SNV.vcf"

## Calculate TiTv

$SAMTOOLS_DIR/bcftools/vcfutils.pl qstats $CORE_PATH/$PROJECT/$SM_TAG/temp/$SM_TAG".UG.EAS.FILTERED.OnBait.TiTv.SNV.vcf" >| \
$CORE_PATH/$PROJECT/$SM_TAG/Reports/TiTv/$SM_TAG"_.titv.txt"

## Eventually, will need a per sample QC report to put into the analysis folder.

rm -rvf $CORE_PATH/$PROJECT/$SM_TAG/temp

echo Sample $8 in Project $1 has finished on `date` >> $CORE_PATH/Sample.Completion.log
