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
Estimate the genome. Write script on R:
```R
data <- read.table("SRR292678.histo")
plot(data[8:600,],type="l")
points(data[1300:1500,])

m <- which(data$V2 == max(data[13:1569,]$V2))

genome_size <- sum(as.numeric(data[2:1569,1]*data[2:1569,2]))/m
# genome_size = 5168217
```
**Genome sie = 5168217**

### Assembling E. coli X genome from paired reads
I have some problems with SPAdes, therefore download [precomputed results](https://disk.yandex.ru/d/4xEI_7gdxzN2D).
Unzip precomputed results.
Than use online version of QUAST to get [contig's and scaffold's reports](https://github.com/rereremin/IB/tree/project3/quast_results)

### Impact of reads with large insert size
Download [precomputed results](https://disk.yandex.ru/d/XHCbTIrvxzN5Y) and unzip like in previous paragraph. 
[Results there](https://github.com/rereremin/IB/tree/project3/quast_results).

Main metrics for single library assembly:
- N50 = 111860
- number of contigs = 210

Main metrics for three-library assembly:
- N50 = 335515
- number of contigs = 105

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
Shiga toxin geens are **stxA**, **stxB**
stxA, stxB are flanked on the 5'-end by 
- Phage DNA adenine methylase
- Phage antitermination protein Q
- Phage protein NinH
- Phage recombination protein NinG
- Gifsy-2 prophage protein
- DNA primase, phage associated
- Putative ATP-dependent helicase

stxA, stxB are flanked on the 3'-end by:
- ORF B78
- Phage head completion protein
- Phage holin/antiholin component S
- Phage lysozyme R
- Phage antirepressor protein
- Phage endopeptidaze Rz

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
bla_1, bla_2 are flanked on the 5'-end by:
- Mobile element protein
bla_1, bla_2 are flanked on the 3'-end by:
- Tryptophan syntase
- Mobile element protein

Also find Class C beta-lactamase gene and Class A beta-lactamase.

Class C beta-lactamase gene is flanked on the 5'-end by:
- Small multidrug resistance (SMR) efflux transporter
Class C beta-lactamase gene is flanked on the 3'-end by:
- Fumarate-reductase group genes

Class A beta-lactamase gene is flanked on the 5'-end by:
- error-prone genes
Class A beta-lactamase gene is flanked on the 3'-end by:
- mobile element protein
  
