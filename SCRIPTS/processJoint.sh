
sample=$1
	
module load vt
module load vcftools
module load bcftools
module load tabix
module load vep/v88
module load zlib/1.2.5


echo "VT STARTED AT: "
echo `date`

### vt normalize
vt normalize $sample.vcf.gz -r /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -o $sample.norm.vcf

### vt decompose 
vt decompose -s $sample.norm.vcf -o $sample.norm.decomp.vcf
bgzip $sample.norm.decomp.vcf
tabix -p vcf $sample.norm.decomp.vcf.gz

### QC clean
vcftools --gzvcf $sample.norm.decomp.vcf.gz --hwe 0.0000001 --min-alleles 2 --max-alleles 2 --minGQ 20 --min-meanDP 10 --max-missing 0.90 --recode --recode-INFO-all --out $sample.norm.decomp
module purge

### Add RSID
module load gatk/4.0.2.1
echo "VariantAnnotator STARTED AT: "
echo `date`
gatk VariantAnnotator \
	-R /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
	-V $sample.norm.decomp.recode.vcf \
	-O $sample.norm.decomp.recode.dbSNP.vcf.gz \
	--dbsnp /rgs01/project_space/cab/Control/common/reference/GATK_support/dbsnp_146.hg38.vcf.gz

module purge

### Fix encoding in VCF prior to ANNOVAR
zcat $sample.norm.decomp.recode.dbSNP.vcf.gz | sed 's/<\*:DEL>/<DEL>/g' > $sample.norm.decomp.recode.dbSNP.sed.vcf
	
module load annovar/20191024

echo "ANNOVAR STARTED AT: "
echo `date`

### Add refGene and gnomad annotations
table_annovar.pl $sample.norm.decomp.recode.dbSNP.sed.vcf /research/rgs01/project_space/cab/Control/common/reference/annovar/humandb -buildver hg38 \
-out $sample.sed -remove \
-nastring . -vcfinput \
-protocol refGene,gnomad211_exome1,gnomad30_genome1 \
-operation g,f,f \
--argument "--splicing_threshold 4, , " 

module purge

echo "SLIVAR STARTED AT: "
echo `date`

### Filter out variants that have MAF >0.05
/research/rgs01/project_space/cab/Control/common/reference/slivar/mne/slivar expr \
	--pass-only \
	--vcf $sample.sed.hg38_multianno.vcf \
	--out-vcf $sample.sed.hg38_multianno.slivar05.vcf \
	--info "variant.call_rate > 0.9 && \
		(INFO.gnomad211_exome_AF=='.' && INFO.gnomad30_genome_AF=='.') || \
		((INFO.gnomad30_genome_AF < 0.05 && \
		INFO.gnomad30_genome_AF_afr < 0.05 && \
		INFO.gnomad30_genome_AF_amr < 0.05 && \
		INFO.gnomad30_genome_AF_fin < 0.05 && \
		INFO.gnomad30_genome_AF_nfe < 0.05) || \
		(INFO.gnomad211_exome_AF < 0.05 && \
		INFO.gnomad211_exome_AF_afr < 0.05 && \
		INFO.gnomad211_exome_AF_amr < 0.05 && \
		INFO.gnomad211_exome_AF_asj < 0.05 && \
		INFO.gnomad211_exome_AF_eas < 0.05 && \
		INFO.gnomad211_exome_AF_fin < 0.05 && \
		INFO.gnomad211_exome_AF_nfe < 0.05 && \
		INFO.gnomad211_exome_AF_sas < 0.05) ) " \

### Compress and Index
module load tabix
bgzip -f $sample.sed.hg38_multianno.slivar05.vcf
tabix -p vcf $sample.sed.hg38_multianno.slivar05.vcf.gz
