# Lab report 
## Differential RNA expression analysis

1. Download input files:
- [SRR941816: fermentation 0 minutes replicate 1](tp.sra.ebi.ac.uk/vol1/fastq/SRR941/SRR941816/SRR941816.fastq.gz)
- [SRR941817: fermentation 0 minutes replicate 2](ftp.sra.ebi.ac.uk/vol1/fastq/SRR941/SRR941817/SRR941817.fastq.gz)
- [SRR941818: fermentation 30 minutes replicate 1](ftp.sra.ebi.ac.uk/vol1/fastq/SRR941/SRR941818/SRR941818.fastq.gz)
- [SRR941819: fermentation 30 minutes replicate 2](ftp.sra.ebi.ac.uk/vol1/fastq/SRR941/SRR941819/SRR941819.fastq.gz)

And files with reference genome and annoteion file from NCBI:
- Saccharomyces cerevisiae strain S288c
- assembly R64
  
2. Aligning with HISAT2
  Unzip files and build genome index:
  ```bash
  hisat2-build GCF_000146045.2_R64_genomic.fna genome_index
  ```
  Run hisat2 in single-end mode (align reads on reference)
  ```bash
  hisat2 -p 16 -x ~/Desktop/BI/learning/bioinf_practise/project6/work/genome_index -U raw_data/SRR941816.fastq | samtools sort > SRR941816.bam
  hisat2 -p 16 -x ~/Desktop/BI/learning/bioinf_practise/project6/work/genome_index -U raw_data/SRR941817.fastq | samtools sort > SRR941817.bam
  hisat2 -p 16 -x ~/Desktop/BI/learning/bioinf_practise/project6/work/genome_index -U raw_data/SRR941818.fastq | samtools sort > SRR941818.bam
  hisat2 -p 16 -x ~/Desktop/BI/learning/bioinf_practise/project6/work/genome_index -U raw_data/SRR941819.fastq | samtools sort > SRR941819.bam
  ```

3. Quantifying with featureCounts
  ```bash
  featureCounts -g gene_id -a annotation.gtf -o SRR SRR941816.bam SRR941817.bam SRR941818.bam SRR941819.bam
  cat SRR | cut -f 1,7-10 > simple_counts_SRR.txt
  ```
  Results of featureCounts in `results/simple_counts_SRR` [there](https://github.com/rereremin/IB/tree/project6/feature_counts)
  
4. Find differentially expressed genes with Deseq2
   Run `deseq2.r` to count metrics:
   ```bash
   cat simple_counts_SRR.txt | R -f scripts/deseq2.r
   ```
   results in files `norm-matrix-deseq.txt` and `result.txt`[there](https://github.com/rereremin/IB/tree/project6/results)

5. Draw heatmap
   ```bash
   cat norm-matrix-deseq2.txt | R -f ../scripts/draw-heatmap.r
   ```
   results files [there](https://github.com/rereremin/IB/tree/project6/results)

6. [link](https://www.yeastgenome.org/goSlimMapper?genes=&uploadFile=genes.txt&slim_type=Yeast+GO-Slim%3A+process&submit=Submit+Form)
