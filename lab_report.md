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
- in SRR292678_R1.fastq.gz **5499346** reads
- in SRR292678_R2.fastq.gz **5499346** reads
- in SRR292770_R1.fastq.gz **5102041** reads
- in SRR292770_R2.fastq.gz **5102041** reads
- in SRR292862_R1.fastq.gz **5102041** reads
- in SRR292862_R2.fastq.gz **5102041** reads

### K-mer profile and genome size estimation

Install jellyfish:
```bash
conda install -c bioconda jellyfish
conda install -c "bioconda/label/cf201901" jellyfish
```
Count k-mers with some flags: -m (length of k-mers), -s (hash with 100 million elements), -t (10 threads), -С (canonical). 
Than make a histogram, which is [here](https://github.com/rereremin/IB/tree/project3/kmers).
```bash
jellyfish count -m 31 -s 100M -t 10 -C <(gzcat SRR292678_R1.fastq.gz) <(gzcat SRR292678_R2.fastq.gz)
jellyfish histo -t 10 mer_counts.jf > SRR292678.histo
```
Estimate the genome size using the formula: `Genome_size = T / N`.

T (Total bases) = 55666752

`N = (M * L) / (L - K + 1)`, depth of coverage

M (kmer peak) = 125

K (kmer size) = 31

L (average read length) = 135.041

**Genome_size = 346400**

### Assembling E. coli X genome from paired reads




