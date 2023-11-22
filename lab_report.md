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

Create new directory `fastqc_results` and run FastQC for these six files:
```bash
fastqc -o fastqc_results SRR292678_R1.fastq.gz SRR292678_R2.fastq.gz
fastqc -o fastqc_results SRR292770_R1.fastq.gz SRR292770_R2.fastq.gz
fastqc -o fastqc_results SRR292862_R1.fastq.gz SRR292862_R2.fastq.gz
```
You can see FastQC results [there](https://github.com/rereremin/IB/tree/project3/fastqc_results).
Number of reads in each file:
- 
