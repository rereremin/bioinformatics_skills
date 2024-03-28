# Lab report
## Dead Man’s Teeth. Introduction to metagenomics analysis
### Amplicon sequencing
We used a script written in the R. You shuold find it [here](https://github.com/rereremin/BI/tree/project7/scripts)

In one of the steps, a table with information on fastq files was obtained. This table contains the total number of reeds in the file, the number of filtered reeds and the number of non-chimeric reeds. Table [here](https://github.com/rereremin/BI/tree/project7/results).

As a result of executing this script we have 4 files, which we need for [MicrobiomeAnalyst](https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/ModuleView.xhtml). Results of this step [here](https://github.com/rereremin/BI/tree/project7/results) 

[Taxa barplot: ](https://github.com/rereremin/BI/blob/project7/results/taxa_alpha_2.svg))

### Shotgun sequencing
1. Shotgun sequence data profiling
   Download report and taxonomic prediction for each sequence.
2. Visualization of the Kraken results as a Sankey diagram
   [Results](https://github.com/rereremin/BI/blob/project7/results/sankey-SRR957742_reads.report.html)
3. Comparison with ancient Tannerella forsythia genome

   Download data for the T. forsythia (genome in fasta and annotation in gff3)

   ```bash
   bwa index GCA_000238215.1_ASM23821v1_genomic.fna
   bwa mem GCA_000238215.1_ASM23821v1_genomic.fna G12_assembly.fna.gz > alignment.sam
   samtools view -S -b alignment.sam > alignment.bam
   samtools sort alignment.bam -o alignment_sorted.bam
   samtools index alignment_sort.bam
   ```
   Then visualized with IGV.
