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
`ext.proteins.fasta` includes 34 unique proteins.

### Localization prediction
Using TargretP we find *other* category of proteins. List with 21 *other* protein:
- g10513.t1
- g10514.t1
- g11513.t1
- g11806.t1
- g11960.t1
- g12510.t1
- g14472.t1
- g15484.t1
- g16318.t1
- g16368.t1
- g2203.t1
- g3428.t1
- g4106.t1
- g4970.t1
- g5237.t1
- g5443.t1
- g5510.t1
- g5927.t1
- g7861.t1
- g8100.t1
- g8312.t1

Firstly make `cur.proteins.sseqid.txt` with list of proteins. Than make fasta with these proteins:
```bash
seqtk subseq augustus.whole.aa cur.proteins.sseqid.txt > cur.ext.proteins.fasta
```

### BLAST search
- g10513.t1 (No significant similarity found)
- g10514.t1 (No significant similarity found)
- **g11513.t1**

  Q96Q05.2, 4e-81,	27.04%, 76%

  RecName: Full=Trafficking protein particle complex subunit 9;

  AltName: Full=NIK- and IKBKB-binding protein;

  AltName: Full=Tularik gene 1 protein [Homo sapiens]

- g11806.t1 (No significant similarity found)
- **g11960.t1**

  Q5DTM8.2, 2e-82, 26.42%, 98%

  RecName: Full=E3 ubiquitin-protein ligase BRE1A; Short=BRE1-A;

  AltName: Full=RING finger protein 20;

  AltName: Full=RING-type E3 ubiquitin transferase BRE1A [Mus musculus]

- g12510.t1 (No significant similarity found)
- **g14472.t1**

  P0DOW4.1, 0.0, 100.00%, 100%

  RecName: Full=Damage suppressor protein [Ramazzottius varieornatus]
  
- **g15484.t1**

  Q155U0.1, 0.0, 45.03%, 78%

  RecName: Full=Vacuolar protein sorting-associated protein 51 homolog;

  AltName: Full=Protein fat-free [Danio rerio]

- g16318.t1 (No significant similarity found)
- g16368.t1 (No significant similarity found)
- **g2203.t1**
  
  Q69ZQ1.2, 2e-126, 35.93%, 75%
  
  RecName: Full=Myogenesis-regulating glycosidase;
  
  AltName: Full=Nuclear envelope transmembrane protein 37;
  
  AltName: Full=Uncharacterized family 31 glucosidase KIAA1161 [Mus musculus]
- **g3428.t1**
  
  P02609.2, 8e-50, 48.21%, 97%

  RecName: Full=Myosin regulatory light chain 11;
  
  AltName: Full=DTNB;

  AltName: Full=Fast skeletal myosin light chain 2; Short=MLC-2;

  AltName: Full=G2;

  AltName: Full=LC2f;

  AltName: Full=Myosin light chain 11;

  AltName: Full=Myosin regulatory light chain 2, skeletal muscle isoform [Gallus gallus]
  
- g4106.t1 (No significant similarity found)
- **g4970.t1**

  Q9Y5Q5.2, 2e-11, 22.96%, 65%

  RecName: Full=Atrial natriuretic peptide-converting enzyme;

  AltName: Full=Corin;

  AltName: Full=Heart-specific serine proteinase ATC2;

  AltName: Full=Pro-ANP-converting enzyme;

  AltName: Full=Transmembrane protease serine 10;

  Contains: RecName: Full=Atrial natriuretic peptide-converting enzyme, N-terminal propeptide;

  Contains: RecName: Full=Atrial natriuretic peptide-converting enzyme, activated protease fragment;

  Contains: RecName: Full=Atrial natriuretic peptide-converting enzyme, 180 kDa soluble fragment;

  Contains: RecName: Full=Atrial natriuretic peptide-converting enzyme, 160 kDa soluble fragment;

  Contains: RecName: Full=Atrial natriuretic peptide-converting enzyme, 100 kDa soluble fragment [Homo sapiens]

- g5237.t1 (No significant similarity found)
- g5443.t1 (No significant similarity found)
- g5510.t1 (No significant similarity found)
- **g5927.t1**
  
  Q9VAI0.1, 3e-17, 33.95%, 17%

  RecName: Full=Probable glucosamine 6-phosphate N-acetyltransferase;

  AltName: Full=Phosphoglucosamine acetylase;

  AltName: Full=Phosphoglucosamine transacetylase [Drosophila melanogaster]

- **g7861.t1**

  B4F769.1, 2e-71, 37.21%, 99%

  RecName: Full=SWI/SNF-related matrix-associated actin-dependent regulator of chromatin subfamily A-like protein 1;

  AltName: Full=HepA-related protein;

  AltName: Full=Sucrose nonfermenting protein 2-like 1 [Rattus norvegicus]

- **g8100.t1**

  Q9VYF2.1, 1e-43,	32.84%, 25%

  RecName: Full=Putative inositol monophosphatase 3; Short=IMP 3; Short=IMPase 3;

  AltName: Full=Inositol-1(or 4)-monophosphatase 3;

  AltName: Full=Myo-inositol monophosphatase A3 [Drosophila melanogaster]
  
- **g8312.t1**

  P38959.2, 2e-35, 22.15%, 93%

  RecName: Full=Vacuolar protein sorting-associated protein 41;

  AltName: Full=Vacuolar morphogenesis protein 2 [Saccharomyces cerevisiae S288C]

  Other blast results [here](https://github.com/rereremin/IB/tree/project4/blast_res)
