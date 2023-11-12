# Project 2

## Creating directories
```python
mkdir project2
cd project2
mkdir raw_data
mkdir working_dir
```
### Downloading raw data
```python
snakemake --cores=all -p roommate.fastq.gz
```
Downloaded reference HA gene from [here](https://www.ncbi.nlm.nih.gov/nuccore/KF848938.1?report=fasta).

### Inspecting roommate data
Count reads in roommate sample:
```python
wc -l roommate.fastq

# 1433060 roommate.fastq
```
Inspect reads with FastQC:
```python
fastqc -o . roommate.fastq
```
### Alignment roommate's data to reference
Create sorted bam file:
```python
snakemake --cores=all -p reference.roommate.sorted.bam
```
Count aligned reads:
```python
samtools flagstat reference.roommate.sorted.bam

# 361349 + 0 in total (QC-passed reads + QC-failed reads)
# 358265 + 0 primary
# 0 + 0 secondary
# 3084 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 361116 + 0 mapped (99.94% : N/A)
# 358032 + 0 primary mapped (99.93% : N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)
```
### Looking for common variants with VarScan

Compute reference length:
```python
cat reference.fasta | grep -v '^>' | wc -m
# 1690
```
Compute number of reads:
```python
awk '{s++}END{print s/4}' roommate.fastq
# 357171
```
Compute average read length:
```python
awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' roommate.fastq
# 145.042
```
Create my.mpileup and use -d flag to increase depth to 60000
```python
samtools mpileup -d 60000 -f reference.fasta reference.roommate.sorted.bam > my.mpileup
```
Run VarScan with N = 0.95
```python
varscan mpileup2snp my.mpileup --min-var-freq 0.95 --variants --output-vcf 1 > VarScan_results_0_95.vcf

# Only SNPs will be reported
# Warning: No p-value threshold provided, so p-values will not be calculated
# Min coverage:	8
# Min reads2:	2
# Min var freq:	0.95
# Min avg qual:	15
# P-value thresh:	0.01
# Reading input from my.mpileup
# 1665 bases in pileup file
# 5 variant positions (5 SNP, 0 indel)
# 0 were failed by the strand-filter
# 5 variant positions reported (5 SNP, 0 indel)
```
5 variants are reported back:
```python
cat VarScan_results_0_95.vcf | awk 'NR>24 {print $1, $2, $3, $4, $5}'

# KF848938.1 72 . A G
# KF848938.1 117 . C T
# KF848938.1 774 . T C
# KF848938.1 999 . C T
# KF848938.1 1260 . A C
```
KF848938.1 72 A -> G, Thr -> Thr, synonymous 

KF848938.1 117 C -> T, Ala -> Ala, synonymous

KF848938.1 774 T -> C, Phe -> Phe, synonymous

KF848938.1 999 C -> T, Gly -> Gly, synonymous

KF848938.1 1260 A -> C, Leu -> Leu, synonymous

These mutations do not contribute to the infection.

### Looking for rare proteins with VarScan

Run VarScan with N = 0.001
```python
varscan mpileup2snp my.mpileup --min-var-freq 0.001 --variants --output-vcf 1 > VarScan_results_0_001.vcf

# Scan_results_0_001.vcf
# Only SNPs will be reported
# Warning: No p-value threshold provided, so p-values will not be calculated
# Min coverage:	8
# Min reads2:	2
# Min var freq:	0.001
# Min avg qual:	15
# P-value thresh:	0.01
# Reading input from my.mpileup
# 1665 bases in pileup file
# 23 variant positions (21 SNP, 2 indel)
# 0 were failed by the strand-filter
# 21 variant positions reported (21 SNP, 0 indel)
```
21 variants are reported back:
```python
cat VarScan_results_0_95.vcf | awk 'NR>24 {print $1, $2, $3, $4, $5}'

# KF848938.1 72 . A G
# KF848938.1 117 . C T
# KF848938.1 254 . A G
# KF848938.1 276 . A G
# KF848938.1 307 . C T
# KF848938.1 340 . T C
# KF848938.1 389 . T C
# KF848938.1 691 . A G
# KF848938.1 722 . A G
# KF848938.1 744 . A G
# KF848938.1 774 . T C
# KF848938.1 802 . A G
# KF848938.1 859 . A G
# KF848938.1 915 . T C
# KF848938.1 999 . C T
# KF848938.1 1043 . A G
# KF848938.1 1086 . A G
# KF848938.1 1213 . A G
# KF848938.1 1260 . A C
# KF848938.1 1280 . T C
# KF848938.1 1458 . T C
```
### Inspecting control data sets

```python
snakemake --cores=all -p control1.fastq
snakemake --cores=all -p control2.fastq
snakemake --cores=all -p control3.fastq
```
```python
wc -l control1.fastq
wc -l control2.fastq
wc -l control3.fastq

# 1026344 control1.fastq
# 933308 control2.fastq
# 999856 control3.fastq
```
Inspect reads with FastQC:
```python
fastqc -o . control1.fastq
fastqc -o . control2.fastq
fastqc -o . control3.fastq
```
Create sorted bam files:
```python
snakemake --cores=all -p reference.control1.sorted.bam
snakemake --cores=all -p reference.control2.sorted.bam
snakemake --cores=all -p reference.control3.sorted.bam
```
Count aligned reads:
```python
samtools flagstat reference.control1.sorted.bam

# 256744 + 0 in total (QC-passed reads + QC-failed reads)
# 256586 + 0 primary
# 0 + 0 secondary
# 158 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 256658 + 0 mapped (99.97% : N/A)
# 256500 + 0 primary mapped (99.97% : N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)
```
```python
samtools flagstat reference.control2.sorted.bam

# 233451 + 0 in total (QC-passed reads + QC-failed reads)
# 233327 + 0 primary
# 0 + 0 secondary
# 124 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 233375 + 0 mapped (99.97% : N/A)
# 233251 + 0 primary mapped (99.97% : N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)
```
```python
samtools flagstat reference.control3.sorted.bam

# 250184 + 0 in total (QC-passed reads + QC-failed reads)
# 249964 + 0 primary
# 0 + 0 secondary
# 220 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 250108 + 0 mapped (99.97% : N/A)
# 249888 + 0 primary mapped (99.97% : N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)
```

### Looking for rare variants in control files

Create my.mpileup files and use -d flag to increase depth to 60000
```python
samtools mpileup -d 60000 -f reference.fasta reference.control1.sorted.bam > my.mpileup1
samtools mpileup -d 60000 -f reference.fasta reference.control2.sorted.bam > my.mpileup2
samtools mpileup -d 60000 -f reference.fasta reference.control3.sorted.bam > my.mpileup3
```

Run VarScan with N = 0.001
```python
varscan mpileup2snp my.mpileup1 --min-var-freq 0.001 --variants --output-vcf 1 > VarScan_results_0_001_1.vcf
# Only SNPs will be reported
# Warning: No p-value threshold provided, so p-values will not be calculated
# Min coverage:	8
# Min reads2:	2
# Min var freq:	0.001
# Min avg qual:	15
# P-value thresh:	0.01
# Reading input from my.mpileup1
# 1665 bases in pileup file
# 58 variant positions (58 SNP, 0 indel)
# 1 were failed by the strand-filter
# 57 variant positions reported (57 SNP, 0 indel)
```
57 variants are reported back:
```python
cat VarScan_results_0_001_1.vcf | awk 'NR>24 {print $1, $2, $3, $4, $5}'
# KF848938.1 38 . T C
# KF848938.1 54 . T C
# KF848938.1 72 . A G
# KF848938.1 95 . A G
# KF848938.1 117 . C T
# KF848938.1 165 . T C
# KF848938.1 183 . A G
# KF848938.1 216 . A G
# KF848938.1 218 . A G
# KF848938.1 222 . T C
# KF848938.1 235 . T C
# KF848938.1 254 . A G
# KF848938.1 276 . A G
# KF848938.1 297 . T C
# KF848938.1 328 . T C
# KF848938.1 340 . T C
# KF848938.1 356 . A G
# KF848938.1 370 . A G
# KF848938.1 389 . T C
# KF848938.1 409 . T C
# KF848938.1 414 . T C
# KF848938.1 421 . A G
# KF848938.1 426 . A G
# KF848938.1 463 . A G
# KF848938.1 516 . A G
# KF848938.1 566 . A G
# KF848938.1 595 . G T
# KF848938.1 597 . A G
# KF848938.1 660 . A G
# KF848938.1 670 . A G
# KF848938.1 691 . A G
# KF848938.1 722 . A G
# KF848938.1 744 . A G
# KF848938.1 774 . T C
# KF848938.1 859 . A G
# KF848938.1 915 . T C
# KF848938.1 987 . A G
# KF848938.1 1008 . T G
# KF848938.1 1031 . A G
# KF848938.1 1043 . A G
# KF848938.1 1056 . T C
# KF848938.1 1086 . A G
# KF848938.1 1089 . A G
# KF848938.1 1213 . A G
# KF848938.1 1260 . A C
# KF848938.1 1264 . T C
# KF848938.1 1280 . T C
# KF848938.1 1281 . T C
# KF848938.1 1286 . T C
# KF848938.1 1339 . T C
# KF848938.1 1358 . A G
# KF848938.1 1398 . T C
# KF848938.1 1421 . A G
# KF848938.1 1460 . A G
# KF848938.1 1482 . A G
# KF848938.1 1580 . T C
# KF848938.1 1591 . T C
```
```python
varscan mpileup2snp my.mpileup2 --min-var-freq 0.001 --variants --output-vcf 1 > VarScan_results_0_001_2.vcf
# Only SNPs will be reported
# Warning: No p-value threshold provided, so p-values will not be calculated
# Min coverage:	8
# Min reads2:	2
# Min var freq:	0.001
# Min avg qual:	15
# P-value thresh:	0.01
# Reading input from my.mpileup2
# 1665 bases in pileup file
# 54 variant positions (54 SNP, 0 indel)
# 2 were failed by the strand-filter
# 52 variant positions reported (52 SNP, 0 indel)
```
52 variants are reported back:
```python
cat VarScan_results_0_001_2.vcf | awk 'NR>24 {print $1, $2, $3, $4, $5}'
# KF848938.1 44 . T C
# KF848938.1 158 . A G
# KF848938.1 165 . T C
# KF848938.1 183 . A G
# KF848938.1 193 . A G
# KF848938.1 216 . A G
# KF848938.1 218 . A G
# KF848938.1 222 . T C
# KF848938.1 254 . A G
# KF848938.1 276 . A G
# KF848938.1 319 . T C
# KF848938.1 340 . T C
# KF848938.1 356 . A G
# KF848938.1 370 . A G
# KF848938.1 398 . A G
# KF848938.1 403 . A G
# KF848938.1 409 . T C
# KF848938.1 414 . T C
# KF848938.1 421 . A G
# KF848938.1 463 . A G
# KF848938.1 499 . A G
# KF848938.1 516 . A G
# KF848938.1 548 . A G
# KF848938.1 591 . A G
# KF848938.1 607 . A G
# KF848938.1 660 . A G
# KF848938.1 670 . A G
# KF848938.1 691 . A G
# KF848938.1 722 . A G
# KF848938.1 744 . A G
# KF848938.1 793 . A G
# KF848938.1 859 . A G
# KF848938.1 898 . A G
# KF848938.1 915 . T C
# KF848938.1 987 . A G
# KF848938.1 1031 . A G
# KF848938.1 1056 . T C
# KF848938.1 1086 . A G
# KF848938.1 1100 . T C
# KF848938.1 1213 . A G
# KF848938.1 1264 . T C
# KF848938.1 1280 . T C
# KF848938.1 1358 . A G
# KF848938.1 1366 . A G
# KF848938.1 1398 . T C
# KF848938.1 1421 . A G
# KF848938.1 1460 . A G
# KF848938.1 1482 . A G
# KF848938.1 1517 . A G
# KF848938.1 1520 . T C
# KF848938.1 1600 . T C
# KF848938.1 1604 . T C
```
```python
varscan mpileup2snp my.mpileup3 --min-var-freq 0.001 --variants --output-vcf 1 > VarScan_results_0_001_3.vcf
# Only SNPs will be reported
# Warning: No p-value threshold provided, so p-values will not be calculated
# Min coverage:	8
# Min reads2:	2
# Min var freq:	0.001
# Min avg qual:	15
# P-value thresh:	0.01
# Reading input from my.mpileup3
# 1665 bases in pileup file
# 61 variant positions (61 SNP, 0 indel)
# 0 were failed by the strand-filter
# 61 variant positions reported (61 SNP, 0 indel)
```
61 variants are reported back:
```python
cat VarScan_results_0_001_3.vcf | awk 'NR>24 {print $1, $2, $3, $4, $5}'
# KF848938.1 38 . T C
# KF848938.1 44 . T C
# KF848938.1 95 . A G
# KF848938.1 105 . A G
# KF848938.1 133 . A G
# KF848938.1 158 . A G
# KF848938.1 165 . T C
# KF848938.1 183 . A G
# KF848938.1 199 . A G
# KF848938.1 216 . A G
# KF848938.1 218 . A G
# KF848938.1 222 . T C
# KF848938.1 228 . T C
# KF848938.1 230 . A G
# KF848938.1 235 . T C
# KF848938.1 254 . A G
# KF848938.1 271 . A G
# KF848938.1 276 . A G
# KF848938.1 297 . T C
# KF848938.1 319 . T C
# KF848938.1 340 . T C
# KF848938.1 356 . A G
# KF848938.1 370 . A G
# KF848938.1 389 . T C
# KF848938.1 409 . T C
# KF848938.1 414 . T C
# KF848938.1 421 . A G
# KF848938.1 463 . A G
# KF848938.1 499 . A G
# KF848938.1 566 . A G
# KF848938.1 597 . A G
# KF848938.1 607 . A G
# KF848938.1 660 . A G
# KF848938.1 670 . A G
# KF848938.1 691 . A G
# KF848938.1 722 . A G
# KF848938.1 744 . A G
# KF848938.1 759 . T C
# KF848938.1 859 . A G
# KF848938.1 915 . T C
# KF848938.1 987 . A G
# KF848938.1 1031 . A G
# KF848938.1 1043 . A G
# KF848938.1 1056 . T C
# KF848938.1 1086 . A G
# KF848938.1 1089 . A G
# KF848938.1 1105 . A G
# KF848938.1 1209 . A G
# KF848938.1 1213 . A G
# KF848938.1 1264 . T C
# KF848938.1 1280 . T C
# KF848938.1 1281 . T C
# KF848938.1 1301 . A G
# KF848938.1 1358 . A G
# KF848938.1 1366 . A G
# KF848938.1 1398 . T C
# KF848938.1 1421 . A G
# KF848938.1 1460 . A G
# KF848938.1 1482 . A G
# KF848938.1 1580 . T C
# KF848938.1 1604 . T C
```
Parse vcf files into txt files:
```python
bcftools query -f '%POS %REF %ALT [ %FREQ]' VarScan_results_0_001.vcf > frequencies.txt
bcftools query -f '%POS %REF %ALT [ %FREQ]' VarScan_results_0_001_1.vcf > frequencies1.txt
bcftools query -f '%POS %REF %ALT [ %FREQ]' VarScan_results_0_001_2.vcf > frequencies2.txt
bcftools query -f '%POS %REF %ALT [ %FREQ]' VarScan_results_0_001_3.vcf > frequencies3.txt
```
Merge these files into one csv file.
Compute mean and standard deviation of frequencies for control samples. 
SNPs with frequencies > 0.5%:

307 C T Pro103Ser - D epitope region
1458 T C Tyr486Tyr
