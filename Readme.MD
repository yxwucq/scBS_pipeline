# snakemake scBS-seq pipeline
Compared with tradtional bs-seq, this pipeline adds a step to remap the unmapped reads to the genome. For detailed explaination, see https://sequencing.qcfail.com/articles/mispriming-in-pbat-libraries-causes-methylation-bias-and-poor-mapping-efficiencies/
The remapping step increased the mapping efficiency about 10%. However, it may introduce some bias about methylation data.
