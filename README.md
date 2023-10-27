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
5. Results in file *VarScan_results_annotated.vcf*
- ID=93043
```
ADP=17;WT=0;HET=0;HOM=1;NC=0;ANN=G|missense_variant|MODERATE|ftsI|b0084|transcript|b0084|protein_coding|1/1|c.1631C>G|p.Ala544Gly|1631/176|

1631/1767|544/588||,G|upstream_gene_variant|MODIFIER|murE|b0085|transcript|b0085|protein_coding||c.-123C>G|||||123|WARNING_TRANSCRIPT_NO_ST

ART_CODON,G|upstream_gene_variant|MODIFIER|murF|b0086|transcript|b0086|protein_coding||c.-1607C>G|||||1607|,G|upstream_gene_variant|MODIFIE

R|mraY|b0087|transcript|b0087|protein_coding||c.-2959C>G|||||2959|,G|upstream_gene_variant|MODIFIER|murD|b0088|transcript|b0088|protein_cod

ing||c.-4044C>G|||||4044|,G|downstream_gene_variant|MODIFIER|cra|b0080|transcript|b0080|protein_coding||c.*4011C>G|||||4011|WARNING_TRANSCR

IPT_NO_START_CODON,G|downstream_gene_variant|MODIFIER|mraZ|b0081|transcript|b0081|protein_coding||c.*2951C>G|||||2951|,G|downstream_gene_va

riant|MODIFIER|rsmH|b0082|transcript|b0082|protein_coding||c.*2008C>G|||||2008|,G|downstream_gene_variant|MODIFIER|ftsL|b0083|transcript|b0

083|protein_coding||c.*1646C>G|||||1646|  GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    1/1:93:18:17:0:17:100%:4,2852E-

10:0:36:0:0:7:10
```
- ID=482698
```
DP=16;WT=0;HET=0;HOM=1;NC=0;ANN=A|missense_variant|MODERATE|acrB|b0462|transcript|b0462|protein_coding|1/1|c.1706A>T|p.Gln569Leu|1706/3150|

1706/3150|569/1049||,A|upstream_gene_variant|MODIFIER|pdeB|b0457|transcript|b0457|protein_coding||c.-4081A>T|||||4081|WARNING_TRANSCRIPT_NO

_START_CODON,A|upstream_gene_variant|MODIFIER|ylaC|b0458|transcript|b0458|protein_coding||c.-3447A>T|||||3447|,A|upstream_gene_variant|MODI

FIER|maa|b0459|transcript|b0459|protein_coding||c.-2780A>T|||||2780|,A|upstream_gene_variant|MODIFIER|hha|b0460|transcript|b0460|protein_co

ding||c.-2390A>T|||||2390|,A|upstream_gene_variant|MODIFIER|tomB|b0461|transcript|b0461|protein_coding||c.-1990A>T|||||1990|,A|upstream_gen

e_variant|MODIFIER|acrR|b0464|transcript|b0464|protein_coding||c.-3063T>A|||||3063|,A|upstream_gene_variant|MODIFIER|mscK|b0465|transcript|

b0465|protein_coding||c.-3838T>A|||||3838|,A|downstream_gene_variant|MODIFIER|acrA|b0463|transcript|b0463|protein_coding||c.*1728A>T|||||17

28|       GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    1/1:87:16:16:0:16:100%:1,6637E-9:0:45:0:0:7:9
```
- ID=852762
```
ADP=14;WT=0;HET=0;HOM=1;NC=0;ANN=G|upstream_gene_variant|MODIFIER|glnH|b0811|transcript|b0811|protein_coding||c.-4758T>C|||||4758|,G|upstre

am_gene_variant|MODIFIER|dps|b0812|transcript|b0812|protein_coding||c.-3851T>C|||||3851|,G|upstream_gene_variant|MODIFIER|rhtA|b0813|transc

ript|b0813|protein_coding||c.-2665T>C|||||2665|,G|upstream_gene_variant|MODIFIER|opgE|b0815|transcript|b0815|protein_coding||c.-165T>C|||||

165|,G|upstream_gene_variant|MODIFIER|mntR|b0817|transcript|b0817|protein_coding||c.-421A>G|||||421|,G|upstream_gene_variant|MODIFIER|ybiR|

b0818|transcript|b0818|protein_coding||c.-885A>G|||||885|,G|upstream_gene_variant|MODIFIER|ybiT|b0820|transcript|b0820|protein_coding||c.-3

201A>G|||||3201|WARNING_TRANSCRIPT_NO_START_CODON,G|downstream_gene_variant|MODIFIER|yliM|b4736|transcript|b4736|protein_coding||c.*2365A>G

|||||2365|,G|downstream_gene_variant|MODIFIER|ompX|b0814|transcript|b0814|protein_coding||c.*1797A>G|||||1797|,G|downstream_gene_variant|MO

DIFIER|mntS|b4705|transcript|b4705|protein_coding||c.*107T>C|||||107|,G|downstream_gene_variant|MODIFIER|ldtB|b0819|transcript|b0819|protei

n_coding||c.*2062T>C|||||2062|,G|intragenic_variant|MODIFIER|rybA|b4416|gene_variant|b4416|||n.852762A>G||||||

GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    1/1:76:14:14:0:14:100%:2,4927E-8:0:36:0:0:8:6
```
- ID=1905761
```
ADP=13;WT=0;HET=0;HOM=1;NC=0;ANN=A|missense_variant|MODERATE|mntP|b1821|transcript|b1821|protein_coding|1/1|c.74G>A|p.Gly25Asp|74/567|74/56

7|25/188||,A|upstream_gene_variant|MODIFIER|yoaE|b1816|transcript|b1816|protein_coding||c.-4176C>T|||||4176|,A|upstream_gene_variant|MODIFI

ER|yoaL|b4751|transcript|b4751|protein_coding||c.-4030C>T|||||4030|,A|upstream_gene_variant|MODIFIER|yobH|b4536|transcript|b4536|protein_co

ding||c.-3164G>A|||||3164|,A|upstream_gene_variant|MODIFIER|yebQ|b1828|transcript|b1828|protein_coding||c.-4515G>A|||||4515|,A|downstream_g

ene_variant|MODIFIER|manX|b1817|transcript|b1817|protein_coding||c.*2742G>A|||||2742|WARNING_TRANSCRIPT_NO_START_CODON,A|downstream_gene_va

riant|MODIFIER|manY|b1818|transcript|b1818|protein_coding||c.*1879G>A|||||1879|,A|downstream_gene_variant|MODIFIER|manZ|b1819|transcript|b1

819|protein_coding||c.*1015G>A|||||1015|,A|downstream_gene_variant|MODIFIER|yobD|b1820|transcript|b1820|protein_coding||c.*502G>A|||||502|,

A|downstream_gene_variant|MODIFIER|rlmA|b1822|transcript|b1822|protein_coding||c.*490C>T|||||490|,A|downstream_gene_variant|MODIFIER|cspC|b

1823|transcript|b1823|protein_coding||c.*1465C>T|||||1465|,A|downstream_gene_variant|MODIFIER|yobF|b1824|transcript|b1824|protein_coding||c

.*1687C>T|||||1687|,A|downstream_gene_variant|MODIFIER|yebO|b1825|transcript|b1825|protein_coding||c.*2500C>T|||||2500|,A|downstream_gene_v

ariant|MODIFIER|mgrB|b1826|transcript|b1826|protein_coding||c.*2862C>T|||||2862|WARNING_TRANSCRIPT_NO_START_CODON,A|downstream_gene_variant

|MODIFIER|kdgR|b1827|transcript|b1827|protein_coding||c.*3547C>T|||||3547|     GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR

1/1:70:13:13:0:13:100%:9,6148E-8:0:44:0:0:11:2
```
- ID=3535147
```
ADP=17;WT=0;HET=0;HOM=1;NC=0;ANN=C|missense_variant|MODERATE|envZ|b3404|transcript|b3404|protein_coding|1/1|c.722T>G|p.Val241Gly|722/1353|7

22/1353|241/450||,C|upstream_gene_variant|MODIFIER|yhgE|b3402|transcript|b3402|protein_coding||c.-2708T>G|||||2708|,C|upstream_gene_variant

|MODIFIER|greB|b3406|transcript|b3406|protein_coding||c.-1665A>C|||||1665|,C|upstream_gene_variant|MODIFIER|yhgF|b3407|transcript|b3407|pro

tein_coding||c.-2238A>C|||||2238|,C|downstream_gene_variant|MODIFIER|hslO|b3401|transcript|b3401|protein_coding||c.*4495A>C|||||4495|,C|dow

nstream_gene_variant|MODIFIER|pck|b3403|transcript|b3403|protein_coding||c.*707A>C|||||707|,C|downstream_gene_variant|MODIFIER|ompR|b3405|t

ranscript|b3405|protein_coding||c.*718T>G|||||718|       GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR

1/1:93:17:17:0:17:100%:4,2852E-10:0:36:0:0:10:7
```
- ID=4390754
```
ADP=15;WT=0;HET=0;HOM=1;NC=0;ANN=T|synonymous_variant|LOW|rsgA|b4161|transcript|b4161|protein_coding|1/1|c.756C>A|p.Ala252Ala|756/1053|756/

1053|252/350||WARNING_TRANSCRIPT_NO_START_CODON,T|upstream_gene_variant|MODIFIER|mscM|b4159|transcript|b4159|protein_coding||c.-1384C>A||||

|1384|WARNING_TRANSCRIPT_NO_START_CODON,T|upstream_gene_variant|MODIFIER|psd|b4160|transcript|b4160|protein_coding||c.-394C>A|||||394|WARNI

NG_TRANSCRIPT_NO_START_CODON,T|upstream_gene_variant|MODIFIER|orn|b4162|transcript|b4162|protein_coding||c.-850G>T|||||850|,T|upstream_gene

_variant|MODIFIER|yjeV|b4670|transcript|b4670|protein_coding||c.-2138G>T|||||2138|,T|upstream_gene_variant|MODIFIER|nnr|b4167|transcript|b4

167|protein_coding||c.-3312G>T|||||3312|,T|upstream_gene_variant|MODIFIER|tsaE|b4168|transcript|b4168|protein_coding||c.-4831G>T|||||4831|,

T|downstream_gene_variant|MODIFIER|yjeO|b4158|transcript|b4158|protein_coding||c.*4736G>T|||||4736|,T|downstream_gene_variant|MODIFIER|queG

|b4166|transcript|b4166|protein_coding||c.*2174C>A|||||2174|       GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR

1/1:81:16:15:0:15:100%:6,4467E-9:0:36:0:0:8:7
```

## Function and role of genes with found SNPs
1. ftsI:
   - a transpeptidase required for synthesis of peptidoglycan in the division septum and is one of several proteins that localize to the septal ring
   - 
2. [acrB](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3933802/)
   - multidrug Efflux Pump
   - belongs to the resistance-nodulation-division
   -  

