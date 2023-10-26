# Laboratory report
## Creating directories
```bash
mkdir project1
mkdir project1/raw_data
 ```

## Downloading files:
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
				
			
		
				
			
		


