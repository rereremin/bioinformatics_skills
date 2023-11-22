# Lab report
## E. coli outbreak investigation, De novo assembly and annotation of bacterial genomes

### Exploring the dataset
Downloading 6 files:
- SRR292678_R1.fastq.gz
- SRR292678_R2.fastq.gz
- SRR292770_R1.fastq.gz
- SRR292770_R2.fastq.gz
- SRR292862_R1.fastq.gz
- SRR292862_R2.fastq.gz

Run FastQC for these files:
```bash
fastqc -o fastqc_results SRR292678_R1.fastq.gz SRR292678_R2.fastq.gz
fastqc -o fastqc_results SRR292770_R1.fastq.gz SRR292770_R2.fastq.gz
fastqc -o fastqc_results SRR292862_R1.fastq.gz SRR292862_R2.fastq.gz
```
You can see FastQC results [there](file:///Users/nikitazherko/Desktop/BI/learning/bioinf_practise/project3/work/fastqc_results/SRR292678_R1_fastqc.html).
