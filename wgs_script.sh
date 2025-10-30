#!/bin/env/bash

set -euo pipefail

# Create conda environment
#conda create -n wgs_malariagen_data -c bioconda

#activate conda
#source /home/waccbip/anaconda3/etc/profile.d/conda.sh
#conda activate wgs_malariagen

CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"

conda activate wgs_malariagen

# Retrieve data_names from ENA repository as TSV files
mkdir -p analysis_files

#-----------------------------------------------(Sequence files)-----------------------------------------------------------------

awk -F'\t' 'NR > 1 {split($7,a,";"); print a[1]; print a[2];}' analysis_files/analysis_report.tsv > download_links.txt

echo "Downloaded ENA metadata to analysis_report.tsv"

cd analysis_files/
# Extracting just first 20  for testing
head -n 20 analysis_report.tsv > twenty_samples.tsv
head -n 10 twenty_samples.tsv

#Separate the links into two files
awk -F'\t' 'NR > 1 {split($7,a,";"); print a[1]; print a[2];}' analysis_report.tsv > download_links.txt

#awk -F'\t' '{split($7,a,";"); print a[1]; print a[2]}' twenty_samples.tsv > download_links.txt

#Create a new directory for bam files and download (Aspear is faster but didnt work)
if [ ! -d bam_files_20 ]; then
   mkdir -p bam_files_20
   cd bam_files_20
fi

while read link; do
    echo "Downloading $link ..."
    wget -c "ftp://${link}"
done < ../download_links.txt


#Confirm files exist at directory
ls -lh bam_files/AR000*.bam
ls -lh bam_files/AR000*.bai

#Confirm files are accurate
samtools quickcheck -v bam_files/*.bam

#Next step is to run piepline using GATK

#Download the Mosquito Reference Genome
wget -O AgamP4.fa.gz https://vectorbase.org/common/downloads/Current_Release/AgambiaePEST/fasta/data/VectorBase-68_AgambiaePEST_Genome.fasta

#Download old version of GATK If you dont hv it, I do so ill proceed
wget https://storage.googleapis.com/gatk-software/gatk-3.7/GenomeAnalysisTK-3.7-0.tar.bz2
tar xjf GenomeAnalysisTK-3.7-0.tar.bz2
export GATK_JAR=$PWD/GenomeAnalysisTK.jar

docker pull broadinstitute/gatk3:3.7-0

#-------------------------------------------------(GATK Analysis Pipeline)-------------------------------------------------------------------------
docker run --rm \
  -v "$PWD/bam_files:/data/bam" \
  -v "$PWD/AgamP4.fa:/data/ref/AgamP4.fa" \
  -v "$PWD/results:/data/out" \
  broadinstitute/gatk3:3.7-0 \
  java -Xmx8g -jar /usr/GenomeAnalysisTK.jar \
  -T RealignerTargetCreator \
  -R /data/ref/AgamP4.fa \
  -I /data/bam/AR0001-C.bam \
  -o /data/out/intervals.list

#create dict files using samtools:
samtools --version


#create samtools index for AgamP4
samtools faidx AgamP4.fa

#create sequence dictionary
docker run --rm \
  -v "$PWD:/data" \
  broadinstitute/gatk3:3.7-0 \
  java -jar /usr/GenomeAnalysisTK.jar \
  -T CreateSequenceDictionary \
  -R /data/AgamP4.fa \
  -o /data/AgamP4.dict

#verify all files needed exist
ls -la AgamP4.fa*


#Run
docker run --rm \
  -v "$PWD:/data/ref" \
  broadinstitute/picard:latest \
  java -jar /usr/picard/picard.jar \
  CreateSequenceDictionary \
  R=/data/ref/AgamP4.fa \
  O=/data/ref/AgamP4.dict

#Modify the BAM file headers to match conig names between your BAM file and the reference genome. The BAM file uses simple contig names like 2R, 3R, 2L while the reference genome uses prefixed names like AgamP4_2L, AgamP4_2R, etc. -no work
docker run --rm \
  -v "$PWD:/data" \
  broadinstitute/picard:latest \
  java -jar /usr/picard/picard.jar \
  ReorderSam \
  I=/data/bam_files/AR0001-C.bam \
  O=/data/bam_files/AR0001-C_reordered.bam \
  SEQUENCE_DICTIONARY=/data/modified_ref/AgamP4_final.dict \
  CREATE_INDEX=true


#Run docker on Indel Aligner using IndelRealigner:
docker run --rm \
  -v "$PWD/bam_files:/data/bam" \
  -v "$PWD/modified_ref:/data/ref" \
  -v "$PWD/results:/data/out" \
  broadinstitute/gatk3:3.7-0 \
  java -Xmx8g -jar /usr/GenomeAnalysisTK.jar \
  -T IndelRealigner \
  -R /data/ref/AgamP4_final.fa \
  -I /data/bam/AR0001-C_reordered.bam \
  -targetIntervals /data/out/intervals.list \
  -o /data/out/AR0001-C_realigned.bam
  
#Run docker on UnifiedGenotyper to generate raw vcf files: 
docker run --rm \
-v "$PWD/bam_files:/data/bam" \
-v "$PWD/modified_ref:/data/ref" \
-v "$PWD/results:/data/out" \
broadinstitute/gatk3:3.7-0 \
java -Xmx8g -jar /usr/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R /data/ref/AgamP4_final.fa \
-I /data/out/AR0001-C_realigned.bam \
-glm BOTH \
-stand_call_conf 30 \
-o /data/out/AR0001-C_raw_variants.vcf


# Check the VCF file was created
ls -lh results/AR0001-C_raw_variants.vcf

# Look at the first few variants
head -20 results/AR0001-C_raw_variants.vcf


#--------------------------------------------------Variant Filtering (Recommended)------------------------------------------------------------
#Filter low-quality variants
docker run --rm \
  -v "$PWD/modified_ref:/data/ref" \
  -v "$PWD/results:/data/out" \
  broadinstitute/gatk3:3.7-0 \
  java -Xmx8g -jar /usr/GenomeAnalysisTK.jar \
  -T VariantFiltration \
  -R /data/ref/AgamP4_final.fa \
  -V /data/out/AR0001-C_raw_variants.vcf \
  --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
  --filterName "GATK_filter" \
  -o /data/out/AR0001-C_filtered_variants.vcf


#-----------------------------------------------Data Analysis in malariagen-----------------------------------------------------------------

pip install malariagen_data scikit-allel pandas matplotlib

#Run wgs_malariagen_script

python ./wgs_malariagen_script


