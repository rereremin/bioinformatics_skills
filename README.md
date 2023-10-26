# Laboratory report
## Creating directories
1) `mkdir project1`
2)  `mkdir project1/raw_data`.

## Downloading files:
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz raw_data
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz
https://figshare.com/articles/dataset/amp_res_2_fastq_zip/10006541/3?file=23769689 download there two files with reads from shotgun sequencing
```
copy file from raw_data into working directory
cp GCF_000005845.2_ASM584v2_genomic.fna.gz ..
cp GCF_000005845.2_ASM584v2_genomic.gff.gz ..
cp amp_res_1.fastq.gz ..
cp amp_res_2.fastq.gz ..
Inspect raw sequencing data with FastQC. Filtering the reads 		

create a virtual environment
conda install -c bioconda fastqc			
CONDA_SUBDIR=osx-64 conda create -n myenv_x86 python=3.9
conda activate myenv_x86
conda config --env --set subdir osx-64


install fastqc
brew install fastqc	
in project1/raw_data/
fastqc -o . amp_res_1.fastq.gz amp_res_2.fastq.gz
in amp_res_1.fastq
Per base sequence quality
Per tile sequence quality\
in amp_res_2.fastq
Per base sequence quality


Filtering the reads 
1) install Trimmomatic
2) Run Trimmomatic in paired end mode, with following parameters
   - `ILLUMINACLIP:TruSeq3-PE.fa:2:30:10` - remove adapters
   - `LEADING:20` - remove leading low quality or N bases (below quality 20)
   - `TRAILING:20` - remove trailing low quality or N bases (below quality 20)
   - `SLIDINGWINDOW:10:20` - scan the read with a 10-base wide sliding window, cutting when the average quality per base drops below 20
   - `MINLEN:20` - drop reads below the 20 bases long
```bash
trimmomatic PE $read1 $read2 trimmed.paired.R1.fastq.gz trimmed.unpaired.R1.fastq.gz trimmed.paired.R2.fastq.gz

trimmed.unpaired.R2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:10:20 MINLEN:20
```

					
				
			
		
				
			
		


