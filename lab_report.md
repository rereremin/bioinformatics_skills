# Lab report
## E. coli outbreak investigation, De novo assembly and annotation of bacterial genomes
Create conda venv.
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

### K-mer profile and genome size estimation (ПЕРЕСЧИТАТЬ)

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
I have some problems with SPAdes, therefore download [precomputed results](https://disk.yandex.ru/d/4xEI_7gdxzN2D).
Unzip precomputed results.
Than use online version of QUAST to get [contig's and scaffold's reports](https://github.com/rereremin/IB/tree/project3/quast_results)

### Impact of reads with large insert size
Download [precomputed results](https://disk.yandex.ru/d/XHCbTIrvxzN5Y) and unzip like in previous paragraph. 
[Results there](https://github.com/rereremin/IB/tree/project3/quast_results).

### Genome Annotation
Download results [here](https://disk.yandex.ru/d/4ZzBdc2IxzZUb) and unzip it. PROKKA results are [here](https://github.com/rereremin/IB/tree/project3/prokka)

### Finding the closest relative of E. coli X
Install Barnapp. Than find rRNA and after that make a `16S.fna` [file](https://github.com/rereremin/IB/tree/project3/prokka). 
```bash 
conda install -c bioconda -c conda-forge barrnap
barrnap --kingdom bac scaffolds.fasta 1> rrna.gff 2> rrna.log
barrnap rrna.gff --outseq rRNA.fna 
```
After that we use “Nucleotide blast” to find a reference genome `Escherichia coli 55989, complete sequence` and save as `55989.fasta` [here](https://github.com/rereremin/IB/tree/project3/prokka).

### Visualisation with Mauve
Add `55989.fasta` and `scaffolds.gbk` and find Shigo toxins type genes.
- stxB (Shiga toxin subunit B precursor)
- stxA (Shiga toxin subunit A precursor)

### Tracing the source of toxin genes 
Shiga toxin geens are stxA, stxB.
These geens are flanked by wrbA, agp, cbpA, cbpM, tor-genes (D, A, C, R, T, S) on the one side and rut-genes (G, F, E, D, C, B_2, A) on the other.

### Antibiotic resistance detection
#### E.coli X
- streptomycin
- amoxicillin
- ampicillin
- cefepime
- cefotaxime
- ceftazidime
- piperacillin
- aztreonam
- ticarcillin
- ceftriaxone
- cephalothin
- sulfamethoxazole
- trimethoprim
- tetracycline
- doxycycline
#### reference genome
- tetracycline
- doxycycline
- minocycline

### Antibiotic resistance mechanism
Genes, which responsible for antiobiotic resistance to β-lactam antibiotics are **bla_1** and **bla_2**.
These genes are flanking by tnpA, luc-genes (A, B, C, D) and lutA on the one side and mobA, mobB, thrB_2, dsbA on the other. 

