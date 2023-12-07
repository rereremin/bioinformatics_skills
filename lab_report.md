# Lab report
## Tardigrades: from genestealers to space marines

### Download files
Download assembled genome GCA_001 949185.1_Rvar_4.0_genomic.fna.gz.

For structural annotation download precomputed AUGUSTUS results - [protein fasta](https://drive.google.com/file/d/1hCEywBlqNzTrIpQsZTVuZk1S9qKzqQAq/view) and [gff](https://drive.google.com/file/d/12ShwrgLkvJIYQV2p1UlXklmxSOOxyxj4/view).

### Structural annotation
We use  `getAnnoFasta.pl` script to extract protein sequencies from gff file. Than count number of proteins:
```bash
perl getAnnoFasta.pl augustus.whole.gff
grep -c "# start gene" augustus.whole.gff
```
Total number of proteins: 16435.

### Physical localization
Using classic BLAST+.

Firstly, install blast in our venv `project4`:
```bash
conda install -c bioconda blast
```
Than make blast-database and perform the search:
```bash
makeblastdb -in augustus.whole.aa -dbtype prot -out blast_database   
blastp -db blast_database -query peptides.fa -outfmt 6 -out blast_res 
```
Make file with extracted proteins:
```bash
cat blast_res | awk '{FS="\t";OFS="\t"} {print $2}' | sort | uniq > peptide.sseqid.txt
# 34 uniq proteins
seqtk subseq augustus.whole.aa peptide.sseqid.txt > ext.proteins.fasta
```

### Localization prediction
