# Discovery of novel predisposing coding and noncoding variants in familial Hodgkin lymphoma
### **Authors:** JE Flerlage, JR Myers, JL Maciaszek, N Oak, SR Rashkin, Y Hui, YD Wang, W Chen, G Wu, TC Chang, K Hamilton, LR Goldin, ML McMaster, M Rotunno, N Caporaso, A Vogt, D Flamish, K Wyatt, J Liu, M Tucker, C Mullighan, KE Nichols, ML Metzger, J Yang, E Rampersaud
### **Abstract:** Familial aggregation of Hodgkin lymphoma (HL) has been demonstrated in large population studies, pointing to genetic predisposition to this hematological malignancy. To understand the genetic variants associated with the development of HL we performed whole genome sequencing on 234 individuals with and without HL from 36 pedigrees that had 2 or more first-degree relatives with HL. A rigorous family-based segregation analysis was performed for identification of coding and noncoding variants using linkage and filtering approaches. Our rule based variant prioritization identified 45 HL risk variants in 28 pedigrees, of which 34 are coding, 11 are noncoding. Overall, there were 4 highly prioritized recurrent variants: a coding variant in KDR (rs56302315), a 5’UTR variant in KLHDC8B (rs387906223), a noncoding variant in an intron of PAX5 (rs147081110), and another noncoding variant in an intron of GATA3 (rs3824666). A newly identified splice variant in KDR (c.3849-2A>C) was observed for one pedigree. Gene-level recurrence of truncating variants was seen for POLR1E in three independent pedigrees as well. While KDR and KLHDC8B have previously been reported, PAX5, GATA3, and POLR1E represent novel observations. There were no known pathogenic variants for eight pedigrees despite high penetrance of HL in the family with our current knowledge of the human genome. Our findings represent a large number of plausible variants based on rigorous analytic pipelines. These findings are only the first step to comprehensively understand the genetic drivers of HL and will grow with time as this dataset is combined with additional pedigrees and as we learn the potential pathogenicity of new variants each day. 

# Analysis Step1:
## Starting from joint VCF with variants for all samples in cohort named 'sjhl.vcf.gz'. The following script was run to pre-process, QC filter, add gnomAD annotations, and remove variants with >5% MAF in any gnomAD sub-population:
```
SCRIPTS/processJoint.sh sjhl
```
# Analysis Step2:
## After creating 1 directory per pedigree with an alias file listing which samples to treat as affected (or obligate carriers) and which samples cannot be carriers (format: #affecteds[TAB]unaffecteds) and pedigree file. The following script was run was run to identify variants that are present in all affected/obligate individuals and not present in unaffected married-in individuals, annotate with additional databases, and convert VCF to tab-delimited format:
```
SCRIPTS/slivar sjhl.sed.hg38_multianno.slivar05.vcf.gz HL*
```
# Analysis Step3:
## After the SLIVAR script finished for all pedigrees we performed coding and noncoding variants. The following script was run to prioritize variants, create IGV snapshots, and pedigree drawings:
```
SCRIPTS/famvp_v2 HL*
```

# Dependencies
## Software:
```
annovar/20191024
bcftools/1.10.2
gatk/4.0.2.1
meme/5.1.0
samtools/1.10
tabix
vcftools
vep/v88
vt
zlib/1.2.5
gcc
R/3.6.1
```
## R libraries:
```
dplyr
readr
stringr
kinship2
```
## ANNOVAR databases:
```
refGene
gnomad211_exome1
gnomad30_genome1
exac03
exac03nontcga
1000g2015aug_all
esp6500siv2_all
dbnsfp35a
intervar_20180118
dbscsnv11
cosmic70
nci60
clinvar_20190305
regsnpintron
regBase
regBase_prediction
wgEncodeRegDnaseClustered
encRegTfbsClustered
wgEncodeBroadHmmGm12878HMM
wgEncodeBroadHmmH1hescHMM
wgEncodeBroadHmmHepg2HMM
wgEncodeBroadHmmHmecHMM
wgEncodeBroadHmmHsmmHMM
wgEncodeBroadHmmHuvecHMM
wgEncodeBroadHmmK562HMM
wgEncodeBroadHmmNhekHMM
wgEncodeBroadHmmNhlfHMM
```

**********

This project is under the general MIT License.

**********
