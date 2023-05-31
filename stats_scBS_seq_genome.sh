#!/bin/bash

sample=$1
ref_length_bed=$2
ref_bed=$3

# shellcheck disable=code

# This code is used to extract the statistics of clean data from json file
# variables: raw_bases, raw_reads, trimmed_bases, trimmed_reads, Q20_bases, Q30_bases
raw_bases=`grep "total_bases" 01.clean_data/${sample}_clean.json | awk 'BEGIN{FS=":|: |,"}{print $2}' | head -n1`
raw_reads=`grep "total_reads" 01.clean_data/${sample}_clean.json | awk 'BEGIN{FS=":|: |,"}{print $2}' | head -n1`
trimmed_bases=`grep "total_bases" 01.clean_data/${sample}_clean.json | awk 'BEGIN{FS=":|: |,"}{print $2}' | head -n2 | tail -n1`
trimmed_reads=`grep "total_reads" 01.clean_data/${sample}_clean.json | awk 'BEGIN{FS=":|: |,"}{print $2}' | head -n2 | tail -n1`
Q20_bases=`grep "q20_bases" 01.clean_data/${sample}_clean.json | awk 'BEGIN{FS=":|: |,"}{print $2}' | head -n1`
Q30_bases=`grep "q30_bases" 01.clean_data/${sample}_clean.json | awk 'BEGIN{FS=":|: |,"}{print $2}' | head -n1`

# Get the mapping efficiency of single-end reads
map_ratio_pe=$(grep "Mapping efficiency:" 02.bismark_map/${sample}_clean_1_bismark_bt2_PE_report.txt | cut -f 2)
map_ratio_se_1=$(grep "Mapping efficiency:" 02.bismark_map/unmap1/${sample}_clean_1.fq.gz_unmapped_reads_1_bismark_bt2_SE_report.txt | cut -f 2)
map_ratio_se_2=$(grep "Mapping efficiency:" 02.bismark_map/unmap2/${sample}_clean_2.fq.gz_unmapped_reads_2_bismark_bt2_SE_report.txt | cut -f 2)

# Get the deduplication ratio of single-end reads
dedup_remove_ratio_pe=$(grep "Total number duplicated alignments removed:" 02.bismark_map/${sample}.deduplication_report.txt | cut -f 2 | sed 's/.*(\([^)]*\)).*/\1/')
dedup_remove_ratio_se_1=$(grep "Total number duplicated alignments removed:" 02.bismark_map/${sample}_R1.deduplication_report.txt | cut -f 2 | sed 's/.*(\([^)]*\)).*/\1/')
dedup_remove_ratio_se_2=$(grep "Total number duplicated alignments removed:" 02.bismark_map/${sample}_R2.deduplication_report.txt | cut -f 2 | sed 's/.*(\([^)]*\)).*/\1/')

# 修改coverage、ratio、depth等函数
# coverage=`bedtools merge -i 02.map/${sample}_mapQ30_sort_rmdup.bed | bedtools intersect -a - -b ${panel_bed} | sort -V -k1,1 -k2,2 -i - | bedtools jaccard -a - -b ${panel_bed} | awk 'BEGIN{FS="\t"}{print $3}' | tail -n1`
# 总的覆盖区域长度
coverage=`bedtools merge -i 04.stat_info/${sample}_sorted_rmdup.bed | bedtools sort -i - | bedtools jaccard -a ${ref_bed} -b - | awk 'BEGIN{FS="\t"}{print $3}' | tail -n1`
# 在Genome上的深度
depth=`bedtools intersect -a ${ref_bed} -b 04.stat_info/${sample}_sorted_rmdup.bed | bedtools genomecov -i - -g ${ref_length_bed} | awk '{FS=OFS="\t"}{ if ($1=="genome"){sum_i+=$2*$3;total_len=$4}}END{print sum_i/total_len}'`

# 甲基化覆盖率
lambda_reads=`grep -c "lambda" 04.stat_info/${sample}_sorted_rmdup.bed`
conversion=`zcat 03.extracted_met_sites/${sample}.CX_report.txt.gz | grep "lambda" | awk 'BEGIN{FS="\t"}{methy+=$4;total+=$4+$5}END{print (1-methy/total)}'`
genome_CG=`zcat 03.extracted_met_sites/${sample}.CX_report.txt.gz | grep "chr" | awk 'BEGIN{FS="\t"}{ if ($6=="CG") {methy+=$4;total+=$4+$5}}END{print methy/total}'`
genome_CHH=`zcat 03.extracted_met_sites/${sample}.CX_report.txt.gz | grep "chr" | awk 'BEGIN{FS="\t"}{ if ($6=="CHH") {methy+=$4;total+=$4+$5}}END{print methy/total}'`
genome_CHG=`zcat 03.extracted_met_sites/${sample}.CX_report.txt.gz | grep "chr" | awk 'BEGIN{FS="\t"}{ if ($6=="CHG") {methy+=$4;total+=$4+$5}}END{print methy/total}'`

echo -e sample'\t'raw_bases'\t'raw_reads'\t'Q20_bases'\t'Q30_bases'\t'trimmed_bases'\t'trimmed_reads'\t'map_ratio_pe'\t'map_ratio_se_1'\t'map_ratio_se_2'\t'dedup_remove_ratio_pe'\t'dedup_remove_ratio_se_1'\t'dedup_remove_ratio_se_2'\t'coverage'\t'depth'\t'lambda_reads'\t'conversion'\t'genome_CG'\t'genome_CHH'\t'genome_CHG > 04.stat_info/${sample}_stats.tsv
echo -e ${sample}'\t'${raw_bases}'\t'${raw_reads}'\t'${Q20_bases}'\t'${Q30_bases}'\t'${trimmed_bases}'\t'${trimmed_reads}'\t'${map_ratio_pe}'\t'${map_ratio_se_1}'\t'${map_ratio_se_2}'\t'${dedup_remove_ratio_pe}'\t'${dedup_remove_ratio_se_1}'\t'${dedup_remove_ratio_se_2}'\t'${coverage}'\t'${depth}'\t'${lambda_reads}'\t'${conversion}'\t'${genome_CG}'\t'${genome_CHH}'\t'${genome_CHG} >> 04.stat_info/${sample}_stats.tsv

echo -e "["$(date)"]\t$sample statistics done!"