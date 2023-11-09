### Questions
**Find out what strains were in this season’s vaccine. Was that one of the flu strains covered by this vaccine**

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
2) Amount of reads:
```bash
awk '{s++}END{print s/4}' roommate.fastq
# 357171
```   
3) Average read length:
```bash
awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' roommate.fastq
# 145.042
```

Create my.mpileup and use -d flag to increase depth to 50000
```bash
samtools mpileup -d 50000 -f reference.fasta roommate_sorted.bam > my.mpileup
```
VarScan
```bash
java -jar VarScan.v2.4.0.jar mpileup2snp my.mpileup --min-var-freq 0.95 --variants --output-vcf 1 > VarScan_results.vcf
```
5 variants are reported back:
```bash
cat VarScan_results.vcf | awk 'NR>24 {print $1, $2, $3, $4, $5}' > SNP.txt
# KF848938.1 72 . A G
# KF848938.1 117 . C T
# KF848938.1 774 . T C
# KF848938.1 999 . C T
# KF848938.1 1260 . A C
```
What do these mutations do?
- KF848938.1 72 A G: T > T (silent), Q > R (missense), N > D (missense)
- KF848938.1 117 C T: A > A (silent), P > L (missense), H > Y (missense)
-  KF848938.1 774 T C: F > F (silent), L > S (missense), STOP > Q (nonsense)
-  KF848938.1 999 C T: G > G (silent), A > V (missense), R > C (missence)
-  KF848938.1 1260 A C: L  > L (silent) , Y > S (missense), M > Q(missemce)
  
**Could they be what allowed your roommate’s virus to escape the antibodies in your body from the flu vaccine?**
**HOW???**

### Look for rare variants with VarScan
Set the minimum variant frequency to 0.001
```bash
java -jar VarScan.v2.4.0.jar mpileup2snp my.mpileup --min-var-freq 0.001 --variants --output-vcf 1 > VarScan_results.vcf
```
Create file with all variants, which reported back (21 SNP).
```bash
cat VarScan_results.vcf | awk 'NR>24 {print $1, $2, $3, $4, $5, $10}' > SNP2.txt

# KF848938.1 72 . A G 1/1:255:16748:16705:6:16698:99,96%:0E0:35:36:4:2:10841:5857
# KF848938.1 117 . C T 1/1:255:20641:20559:36:20523:99,82%:0E0:35:37:27:9:13352:7171
# KF848938.1 276 . A G 0/1:23:36445:36407:36335:63:0,17%:4,2944E-3:37:35:22237:14098:30:33
# KF848938.1 307 . C T 0/1:255:37025:36919:36567:346:0,94%:2,5986E-65:36:35:22110:14457:180:166
# KF848938.1 340 . T C 0/1:21:37344:37188:37120:62:0,17%:7,6822E-3:37:36:22942:14178:38:24
# KF848938.1 389 . T C 0/1:38:31716:31542:31470:68:0,22%:1,2628E-4:37:37:15911:15559:43:25
# KF848938.1 691 . A G 0/1:25:39006:38858:38784:67:0,17%:2,9912E-3:37:34:20977:17807:23:44
# KF848938.1 722 . A G 0/1:38:37614:37570:37489:76:0,2%:1,5322E-4:37:36:20700:16789:39:37
# KF848938.1 744 . A G 0/1:24:37988:37936:37863:65:0,17%:3,5838E-3:37:33:20555:17308:35:30
# KF848938.1 774 . T C 1/1:255:37964:37792:7:37781:99,97%:0E0:32:38:7:0:19560:18221
# KF848938.1 802 . A G 0/1:59:43606:43498:43395:100:0,23%:1,0306E-6:37:35:20200:23195:32:68
# KF848938.1 859 . A G 0/1:24:35543:35390:35326:62:0,18%:3,9585E-3:37:37:14363:20963:27:35
# KF848938.1 915 . T C 0/1:28:36629:36556:36482:67:0,18%:1,4547E-3:35:35:17836:18646:37:30
# KF848938.1 999 . C T 1/1:255:32106:31682:37:31642:99,87%:0E0:36:36:20:17:16936:14706
# KF848938.1 1043 . A G 0/1:24:31188:31147:31088:57:0,18%:3,6512E-3:36:33:17051:14037:19:38
# KF848938.1 1086 . A G 0/1:32:24946:24946:24892:53:0,21%:6,2624E-4:36:35:12474:12418:21:32
# KF848938.1 1213 . A G 0/1:36:24884:24813:24755:56:0,23%:2,2342E-4:37:36:8797:15958:24:32
# KF848938.1 1260 . A C 1/1:255:22785:22754:2:22740:99,94%:0E0:32:37:0:2:9694:13046
# KF848938.1 1280 . T C 0/1:20:23206:23187:23143:43:0,19%:9,2871E-3:37:35:10997:12146:24:19
# KF848938.1 1458 . T C 0/1:255:26811:26710:26486:221:0,83%:3,6055E-40:37:35:6795:19691:79:142
# KF848938.1 1460 . A G 0/1:20:26715:26695:26643:47:0,18%:9,269E-3:36:35:6901:19742:9:38
```

###  Inspect and align the control sample sequencing data
Download control files:
```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/008/SRR1705858/SRR1705858.fastq.gz raw_data
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/009/SRR1705859/SRR1705859.fastq.gz raw_data
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/009/SRR1705859/SRR1705860.fastq.gz raw_data
```
Unzip and rename files into control1.fastq, control2.fastq, control3.fastq

Count reads in control1.fastq, control2.fastq, control3.fastq:
```bash
awk '{s++}END{print s/4}' control1.fastq
# 256586

awk '{s++}END{print s/4}' control2.fastq
# 233327

awk '{s++}END{print s/4}' control3.fastq
# 249964
```
We need calculate average read length in each file and after calculate rough coverage estimate using formula:

*coverage = (read count * read length ) / total genome size*
```bash
awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' control1.fastq
# 148.561

awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' control2.fastq
# 148.446

awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' control3.fastq
# 148.703
```
```
coverage_1 = 256586 * 148.561 / 1690 = 22555

coverage_2 = 233327 * 148.446 / 1690 = 20494

coverage_3 = 249964 * 148.703 / 1690 = 21994
```
Align each control fastq file to the reference:
```bash
bwa mem reference.fasta control1.fastq | samtools view -S -b - | samtools sort -o - > control1_sorted.bam
bwa mem reference.fasta control2.fastq | samtools view -S -b - | samtools sort -o - > control2_sorted.bam
bwa mem reference.fasta control3.fastq | samtools view -S -b - | samtools sort -o - > control3_sorted.bam
```

Use samtools to index each of the new control alignment (bam) files
```bash
samtools index control1_sorted.bam
samtools index control2_sorted.bam
samtools index control3_sorted.bam
```
###  Use VarScan to look for rare variants in the reference files.
Use mpileup to create files 
```bash
samtools mpileup -d 50000 -f reference.fasta control1_sorted.bam > my.mpileup1
samtools mpileup -d 50000 -f reference.fasta control2_sorted.bam > my.mpileup2
samtools mpileup -d 50000 -f reference.fasta control3_sorted.bam > my.mpileup3
```
Create vcf files
```bash
java -jar VarScan.v2.4.0.jar mpileup2snp my.mpileup1 --min-var-freq 0.001 --variants --output-vcf 1 > VarScan_results_1.vcf
java -jar VarScan.v2.4.0.jar mpileup2snp my.mpileup2 --min-var-freq 0.001 --variants --output-vcf 1 > VarScan_results_2.vcf
java -jar VarScan.v2.4.0.jar mpileup2snp my.mpileup3 --min-var-freq 0.001 --variants --output-vcf 1 > VarScan_results_3.vcf
```
Parse vcf files:
```bash
cat VarScan_results_1.vcf | awk 'NR>24 {print $1, $2, $3, $4, $5, $10}' > VS1.txt
cat VarScan_results_2.vcf | awk 'NR>24 {print $1, $2, $3, $4, $5, $10}' > VS2.txt
cat VarScan_results_3.vcf | awk 'NR>24 {print $1, $2, $3, $4, $5, $10}' > VS3.txt
```
