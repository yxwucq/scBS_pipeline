import os 

file_list = os.listdir("00.raw_data")
name_list = set([x.split('.')[0][:-2] for x in file_list])
gatk_ref_path = "/datb1/wuyuxuan/resources/GATK4_hg38" # includes GATK4_hg38 data
bismark_ref_path = "/datb1/wuyuxuan/resources/bismark_hg38_with_lambda" # includes bismark data
merge_bismark_cov_path = "/home/wuyuxuan/py_scripts/merge_bismark_cov.py"

"""
onsuccess:
    shell("python /home/wuyuxuan/py_scripts/notify.py SampleMetDone! sample done")
onerror:
    shell("python /home/wuyuxuan/py_scripts/notify.py Error! an error occured")
"""

rule all:
    input:
        CX_report=expand("03.extracted_met_sites/{sample}.CX_report.txt.gz", sample=name_list),
        stats=expand("04.stat_info/{sample}_met_stats.tsv", sample=name_list)
        

## step1 - filter
rule TrimData:
    input:
        R1="00.raw_data/{sample}_1.fq.gz",
        R2="00.raw_data/{sample}_2.fq.gz"
    output:
        pre_clean_R1="01.clean_data/{sample}_pre_clean_1.fq.gz",
        pre_clean_R2="01.clean_data/{sample}_pre_clean_2.fq.gz",
        clean_R1="01.clean_data/{sample}_clean_1.fq.gz",
        clean_R2="01.clean_data/{sample}_clean_2.fq.gz",
        report_json="01.clean_data/{sample}_clean.json",
        report_html="01.clean_data/{sample}_clean.html"
    threads: 4
    conda:
        "tgspipe"
    shell: """
        cutadapt -u 6 {input.R1} | gzip -c > {output.pre_clean_R1}
        cutadapt -u 19 {input.R2} | gzip -c > {output.pre_clean_R2}
        fastp -D -i {output.pre_clean_R1} -I {output.pre_clean_R2} -o {output.clean_R1} -O {output.clean_R2} \
        -j {output.report_json} -h {output.report_html} 
        echo -e "["$(date)"]\\t{wildcards.sample} QC done!"
    """

## step2 - bismark map
rule BS_mapping:
    input:
        clean_R1="01.clean_data/{sample}_clean_1.fq.gz",
        clean_R2="01.clean_data/{sample}_clean_2.fq.gz",
    output: 
        # bam file
        bis_rmdup_bam="02.bismark_map/{sample}.deduplicated.bam",
        bis_rmdup_bam_R1="02.bismark_map/{sample}_R1.deduplicated.bam",
        bis_rmdup_bam_R2="02.bismark_map/{sample}_R2.deduplicated.bam",
        # report the duplication rate
        bis_rmdup_report="02.bismark_map/{sample}.deduplication_report.txt",
        bis_rmdup_report_1="02.bismark_map/{sample}_R1.deduplication_report.txt",
        bis_rmdup_report_2="02.bismark_map/{sample}_R2.deduplication_report.txt",
        # report the mapping rate
        bis_mapping_report="02.bismark_map/{sample}_clean_1_bismark_bt2_PE_report.txt",
        bis_mapping_report_1="02.bismark_map/unmap1/{sample}_clean_1.fq.gz_unmapped_reads_1_bismark_bt2_SE_report.txt",
        bis_mapping_report_2="02.bismark_map/unmap2/{sample}_clean_2.fq.gz_unmapped_reads_2_bismark_bt2_SE_report.txt",
    threads: 4
    conda:
        "bs_map"
    shadow: "full"
    shell: """
        bismark -p {threads} --parallel 2 --fastq --pbat --unmapped --output_dir ./02.bismark_map {bismark_ref_path} -1 {input.clean_R1} -2 {input.clean_R2}
        bismark -p {threads} --parallel 2 --fastq --pbat --unmapped --output_dir ./02.bismark_map/unmap1 {bismark_ref_path} ./02.bismark_map/{wildcards.sample}_clean_1.fq.gz_unmapped_reads_1.fq.gz
        bismark -p {threads} --parallel 2 --fastq --unmapped --output_dir ./02.bismark_map/unmap2 {bismark_ref_path} ./02.bismark_map/{wildcards.sample}_clean_2.fq.gz_unmapped_reads_2.fq.gz

        deduplicate_bismark -p --output_dir ./02.bismark_map -o {wildcards.sample} --bam 02.bismark_map/{wildcards.sample}_clean_1_bismark_bt2_pe.bam
        deduplicate_bismark -s --output_dir ./02.bismark_map -o {wildcards.sample}_R1 --bam ./02.bismark_map/unmap1/{wildcards.sample}_clean_1.fq.gz_unmapped_reads_1_bismark_bt2.bam
        deduplicate_bismark -s --output_dir ./02.bismark_map -o {wildcards.sample}_R2 --bam ./02.bismark_map/unmap2/{wildcards.sample}_clean_2.fq.gz_unmapped_reads_2_bismark_bt2.bam
    """

## step3 - extract methylation
rule ExtractMethyl:
    input: 
        bis_rmdup_bam="02.bismark_map/{sample}.deduplicated.bam",
        bis_rmdup_bam_R1="02.bismark_map/{sample}_R1.deduplicated.bam",
        bis_rmdup_bam_R2="02.bismark_map/{sample}_R2.deduplicated.bam",
    output:
        merged_coverage="03.extracted_met_sites/{sample}_merge.cov.gz",
        CX_report="03.extracted_met_sites/{sample}.CX_report.txt.gz",
    threads: 4
    conda:
        "bs_map"
    shadow: "full"
    shell: """
    bismark_methylation_extractor --bedgraph --gzip --multicore {threads} --o ./03.extracted_met_sites --buffer_size 8G --genome_folder {bismark_ref_path} {input.bis_rmdup_bam}
    bismark_methylation_extractor -s --gzip --bedgraph --multicore {threads} --o ./03.extracted_met_sites --buffer_size 8G --genome_folder {bismark_ref_path} {input.bis_rmdup_bam_R1} {input.bis_rmdup_bam_R2}

    python {merge_bismark_cov_path} -o {wildcards.sample} --dir 03.extracted_met_sites \
    03.extracted_met_sites/{wildcards.sample}.deduplicated.bismark.cov.gz 03.extracted_met_sites/{wildcards.sample}_R1.deduplicated.bismark.cov.gz 03.extracted_met_sites/{wildcards.sample}_R2.deduplicated.bismark.cov.gz
    coverage2cytosine --dir 03.extracted_met_sites --gzip --CX --genome_folder {bismark_ref_path} -o {wildcards.sample} {output.merged_coverage}
    """

## step4 - stats info
rule Stats:
    input:
        # raw reads
        report_json="01.clean_data/{sample}_clean.json",
        report_html="01.clean_data/{sample}_clean.html",
        # mapping rate
        bis_rmdup_bam="02.bismark_map/{sample}.deduplicated.bam",
        bis_rmdup_bam_R1="02.bismark_map/{sample}_R1.deduplicated.bam",
        bis_rmdup_bam_R2="02.bismark_map/{sample}_R2.deduplicated.bam",
        # duplication rate
        bis_rmdup_report="02.bismark_map/{sample}.deduplication_report.txt",
        # methylation rate
        CX_report="03.extracted_met_sites/{sample}.CX_report.txt.gz",
        # genome coverage
        chrom_sizes="ref/hg38.chrom.sizes",
        ref_bed="ref/hg38.genome.bed",
        script="stats_scBS_seq_genome.sh",
    output:
        stats="04.stat_info/{sample}_met_stats.tsv",
        sorted_rmdup_bam="04.stat_info/{sample}_sorted_rmdup.bam",
        sorted_rmdup_bed="04.stat_info/{sample}_sorted_rmdup.bed",
    threads: 4
    shadow: "full"
    conda:
        "ngspipe"
    shell: """
        samtools merge -@ {threads} 04.stat_info/{wildcards.sample}_merged.bam {input.bis_rmdup_bam} {input.bis_rmdup_bam_R1} {input.bis_rmdup_bam_R2}
        samtools sort -@ {threads} 04.stat_info/{wildcards.sample}_merged.bam -o 04.stat_info/{wildcards.sample}_merged_sorted.bam 
        samtools rmdup 04.stat_info/{wildcards.sample}_merged_sorted.bam 04.stat_info/{wildcards.sample}_merged_rmdup.bam
        samtools sort -@ {threads} 04.stat_info/{wildcards.sample}_merged_rmdup.bam -o {output.sorted_rmdup_bam}
        samtools index -@ {threads} {output.sorted_rmdup_bam}
        bedtools bamtobed -i {output.sorted_rmdup_bam} > {output.sorted_rmdup_bed}
        sh {input.script} {wildcards.sample} {input.chrom_sizes} {input.ref_bed}
    """



