# h+ progect lab journal

Download raw data, `23_and_me.txt` and `genetek.vcf`.

Preprocessed data with plink tool:

```
plink --23file 23_and_me_input.txt --recode vcf --out snps_clean --output-chr MT --snps-only just-acgt
```

### 1. Origins and haplogroups

Identification of maternal mtDNA haplogroup. Uploaded raw 23_and_me.txt file to the 
[James Lick's online tool](https://dna.jameslick.com/mthap/) with default settings, the resulting best matching 
haplogroups are as follows: H(T152C), H1(T152C), H, H16(T152C), H3(T152C), H46, H52, H69, H9. 
The results are saved as html page in the working_directory. 

Same raw 23_and_me.txt file was loaded to [MorleyDNA.com Y-SNP Subclade Predictor](https://isogg.org/wiki/Y-DNA_tools) 
with default settings:

```
These are the most likely classifications for your Y-SNP data:

* R1a1a [R1a-L168 (R1a-M17, R1a-M198)]
* M3 [M-P117 (M-P118)]
* K [K-PF5504 (K-PF5493, K-PF5480)]
* R1b1a2a1a2a1b3~2 [R1b-L421 (R1b-L433, R1b-L88)]
* N1a [N-M96 (N-CTS7095, N-P189)]
```
Detailed results are in the working_dir.

### 2. Annotation - sex and eye color. 

Male sex, after little investigation of 23_and_me.txt file and looking at chromosomes.

Eye color identification, using clean snps file from 23&Me:

step 0: 

`15	28365618	rs12913832	A	G	.	.	PR	GT	0/1` -> A/G -> not blue

step 1: 

`6	396321	rs12203592	C	T	.	.	PR	GT	0/1` -> C/T

`5	33951693	rs16891982	C	G	.	.	PR	GT	0/1` -> C/G

No rs6119471 was found.

Try genotek.vcf file:

`chr15	28365618	rs12913832	A	G	.	.	.	GT	0/1` -> A/G -> not blue

`chr6	396321	rs12203592	C	T	.	.	.	GT	0/1` -> C/T 

`chr5	33951693	rs16891982	C	G	.	.	.	GT	0/1` -> C/G

No rs6119471 was found too.

Switch to step 2:

rs12896399 was found only in 23&Me data:

`14	92773663	rs12896399	G	.	.	.	PR	GT	0/0` -> G/G

Based on these results we fail to predict eye color based on [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3694299/).

Try with [The IrisPlex System](https://hirisplex.erasmusmc.nl).

input: rs12913832 0, rs1800407 2, rs12896399 0, rs16891982 1, rs1393350 0, rs12203592 1. 

output: blue eye, p-value 0.905; intermediate eye p-value 0.083; brown eye p-value 0.012 -> more likely brown

Saved in working_dir as Result.csv.

### 3. Annotation of all SNPs, selection of clinically relevant ones

##### 3.1 SnpEff/SnpSift

`java -jar snpEff.jar GRCh37.75 snps_clean.vcf  > snps_snpeff.vcf`

Download [ClinVar database](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz). 

`java -jar SnpSift.jar annotate clinvar.vcf  snps_clean.vcf > snps_clean_snpsift_clinvar.vcf`

Download [GWAS catalog](https://www.ebi.ac.uk/gwas/api/search/downloads/full). 

`java -jar SnpSift.jar gwasCat -db gwas_catalog.tsv snps_clean.vcf > snps_clean_gwascat.vcf`


##### 3.2 VEP

`awk '($32!="-")' vep_output.txt | grep risk_factor | cut -f 1-3 | sort | uniq`

i3000469	2:138759649-138759649	T

i6007787	2:234183368-234183368	G

i6058143	1:161479745-161479745	G

i6059141	8:133909974-133909974	G

rs1049296	3:133494354-133494354	T - C1 and C2 transferrin ?

rs1169288	12:121416650-121416650	C

rs12150220	17:5485367-5485367	T - slightly increased risk for several autoimmune diseases

rs13266634	8:118184783-118184783	T - increased risk for type-2 diabetes

rs1801197	7:93055753-93055753	G

rs1801274	1:161479745-161479745	G - cancer progression risk

rs1801275	16:27374400-27374400	G

rs1801394	5:7870973-7870973	G  - higher risk for meningiomas

rs1801968	9:132580901-132580901	G

rs2241880	2:234183368-234183368	G - increased risk for Crohn's disease in Caucasians

rs231775	2:204732714-204732714	G - increased risk for thyroiditis

rs460897	1:196716319-196716319	T

rs4961	4:2906707-2906707	T - 1.8x increased risk for high blood pressure

rs5174	1:53712727-53712727	T - 1.3x increased risk of heart disease

rs61747071	16:53720436-53720436	T

rs6265	11:27679916-27679916	T - ADHA, depression resistance

rs6280	3:113890815-113890815	T

rs6504649	17:48437456-48437456	G

rs699	1:230845794-230845794	G - increased risk of hypertension

More details about found SNPs:

rs1049296	3:133494354-133494354	T - C/T heterozygote.
Heterozygote carrying both C1 and C2 transferrin subtypes; very slightly higher risk for Alzheimers



