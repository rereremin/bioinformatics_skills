# Laboratory report
## Creating directories
```bash
mkdir project1
mkdir project1/raw_data
 ```

## Download files:
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz raw_data
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz raw_data
```
Download [there](https://figshare.com/articles/dataset/amp_res_2_fastq_zip/10006541/3?file=23769689 ) two files with reads from shotgun sequencing. Also using wget.

Copy file from raw_data into working directory
```bash
cp GCF_000005845.2_ASM584v2_genomic.fna.gz ..
cp GCF_000005845.2_ASM584v2_genomic.gff.gz ..
cp amp_res_1.fastq.gz ..
cp amp_res_2.fastq.gz ..
```
Inspect raw sequencing data with FastQC. Filtering the reads 		

## Create virtual environment `project1`
```bash
conda install -c bioconda fastqc			
CONDA_SUBDIR=osx-64 conda create -n project1 python=3.9
conda activate project1
conda config --env --set subdir osx-64
```
Install **FastQC**
```bash
install fastqc
brew install fastqc
```
Run **FastQC** in project1/raw_data/
```bash
fastqc -o . amp_res_1.fastq.gz amp_res_2.fastq.gz
```
Check output html-files:
1) Very unussual results in amp_res_1.fastq: Per base sequence quality, Per tile sequence quality.
2) very unussual results in amp_res_1.fastq: Per base sequence quality.

## Filtering the reads 
1) install Trimmomatic
```bash
conda install -c bioconda trimmomatic
```
   
2) Run Trimmomatic in paired end mode (*PE*), with following parameters
   - `ILLUMINACLIP:TruSeq3-PE.fa:2:30:10` - remove adapters
   - `LEADING:20` - remove leading low quality or N bases (below quality 20)
   - `TRAILING:20` - remove trailing low quality or N bases (below quality 20)
   - `SLIDINGWINDOW:10:20` - scan the read with a 10-base wide sliding window, cutting when the average quality per base drops below 20
   - `MINLEN:20` - drop reads below the 20 bases long
   - use Phred+33
```bash
trimmomatic PE $read1 $read2 trimmed.paired.R1.fastq.gz trimmed.unpaired.R1.fastq.gz trimmed.paired.R2.fastq.gz

trimmed.unpaired.R2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:10:20 MINLEN:20
```
3) Check amount of lines for each output file, devide 4 every result.
```bash
wc -l trimmed.paired.R1.fastq.gz 
wc -l trimmed.paired.R2.fastq.gz
wc -l trimmed.unpaired.R1.fastq.gz
wc -l trimmed.unpaired.R2.fastq.gz
```
Get 
```
37733
37474
784
25
```		
## Aligning sequences to reference
1. Install bwa
```bash
pip install bwa
```
2. Index reference file (*.fna*)
```bash
bwa index GCF_000005845.2_ASM584v2_genomic.fna.gz
```

4. Align reads to reference, use trimmed reads
```bash
bwa mem GCF_000005845.2_ASM584v2_genomic.fna.gz trimmed.paired.R1.fastq.gz trimmed.paired.R2.fastq.gz > alignment.sam
```

5. Compress and sort SAM file
```bash
samtools view -S -b alignment.sam > alignment.bam
```
```bash
samtools flagstat alignment.bam
# 891649 + 0 mapped (99.87% : N/A)
```

6. Sort and index bam file, samtools version = 1.18 (latest)
   
```baash
samtools sort alignment.bam -o alignment.sorted.bam
samtools index alignment.sorted.bam
```

7. I'm working in IGV on my local computeer (MacOS).
Convert reference genome in fasta format

```bash
gzip -dc GCF_000005845.2_ASM584v2_genomic.fna.gz > GCF_000005845.2_ASM584v2_genomic.fnagzip -dc GCF_000005845.2_ASM584v2_genomic.fna.gz > GCF_000005845.2_ASM584v2_genomic.fna
``` 

## Variant calling
1. Create my.mpileup
```bash
samtools mpileup -f GCF_000005845.2_ASM584v2_genomic.fna alignment.sorted.bam > my.mpileup
```

2. We need use VarScan, therefore install it using conda. 
```bash
conda install -c bioconda varscan
```
3. Move VarScan.v2.4.0.jar in our working directory.

4. Create vcf-file with N = 0.2
```bash
java -jar VarScan.v2.4.0.jar mpileup2snp my.mpileup --min-var-freq 0.2 --variants --output-vcf 1 > VarScan_results.vcf
```

## Variant effect prediction
Add two tracks (annotation in gff and vcf file) in IGV. Then find SNPs on the largest scale.
- in position 93043 missense mutation (A -> G) in ftsI
- in position 482698 missense mutation (Q -> L) in acrB
- in position 852762 no genes
- in position 1905761 missense mutation (G -> D) in mntP
- in position 3535147 missense mutation (V -> A) in envZ
- in position 4390754 missense mutation (A -> S) in rsgA


## Automatic SNP annotation
1. We need a SnpEff, therefore install it using conda
   ```bash
   conda install -c bioconda snpeff
   ```
2. Download *GCF_000005845.2_ASM584v2_genomic.gbff.gz* and move it i working directory.
3. Create *snpEff.config* with one line 'k12.genome : ecoli_K12'. Then
   ```bash
   mkdir -p data/k12
   gunzip GCF_000005845.2_ASM584v2_genomic.gbff.gz
   cp GCF_000005845.2_ASM584v2_genomic.gbff data/k12/genes.gbk
   snpEff build -genbank -v k12
   snpEff ann k12 VarScan_results.vcf > VarScan_results_annotated.vcf
   ```
5. Results in file [VarScan_results_annotated.vcf](https://github.com/rereremin/IB/blob/project1/VarScan_results_annotated.vcf)

## Functions and roles of genes with found SNPs
1. ftsI:
   - a transpeptidase required for synthesis of peptidoglycan in the division septum and is one of several proteins that localize to the septal ring
   - check more info [there](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC305773/#:~:text=FtsI%20(also%20called%20PBP3)%20of,localize%20to%20the%20septal%20ring.)
2. acrB
   - multidrug Efflux Pump
   - belongs to the resistance-nodulation-division
   - check more info [there](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3933802/)

