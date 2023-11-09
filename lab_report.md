### Questions
*Find out what strains were in this season’s vaccine. Was that one of the flu strains covered by this vaccine*
### Creating virtual environment
```bash			
CONDA_SUBDIR=osx-64 conda create -n project2 python=3.9
conda activate project2
```
### Creating directories
- raw_data
- results
### Download data
- `wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/001/SRR1705851/SRR1705851.fastq.gz raw_data ` - data from room
- `sequence.fasta` reference [sequence](https://www.ncbi.nlm.nih.gov/nuccore/KF848938.1?report=fasta)
### Before start
```bash
fastqc -o results SRR1705851.fastq.gz 
```
trimmomatic SE
```bash
trimmomatic SE SRR1705851.fastq.gz roommate.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:10:20 MINLEN:20
```
Index reference file 
```bash
bwa index reference.fasta
```
Pipeline for alignment
```bash
bwa mem reference.fasta roommate.fastq| samtools view -S -b - | samtools sort -o - > roommate_sorted.bam
```
Check percentage of mapped reads
```bash
samtools flagstat roommate_sorted.bam 
# 359456 + 0 mapped (99.80% : N/A)
```
### Look common variants with VarScan
Index bam file
```bash
samtools index roommate_sorted.bam
```
#### Count length reference, amount of reads and average read length:
1) Length reference: 
```bash
cat reference.fasta | grep -v '^>' | wc -m
# 1690
```
2) Amount of reads: `360193`
3) Average read length:
```bash
awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' roommate.fastq
# 145.042
```

Create my.mpileup and use -d flag to increase depth to 50000
```bash
samtools mpileup -d 50000 -f reference.fasta roommate_sorted.bam > my.mpileup
```
