#!/bin/bash
#
#

if [ $# -eq 0 ]; then
	echo "Usage: $0 [VariantLIST] [OUTPUT] [VCF] [REFERENCE] [MOTIFDB] [TF_list] [OUTPUT]";
	echo "";
	echo "[VariantLIST] = file containing one variant per line; CHROM<TAB>POS<TAB>ID<TAB>QUAL<TAB>REF<TAB>ALT; subset of the VCF file";
	echo "";
	echo "[OUTPUT] = Output file name. Default = output_FIMOsnp.txt"
	echo "";
	echo "[VCF] = VCF file with variants to generate consensus sequences"
	echo "";
	echo "[REFERENCE] = FASTA of the reference genome";
	echo "";
	echo "[MOTIFDB] = A motif database in meme format - /hpcf/authorized_apps/rhel7_apps/meme/install/5.1.0/db/motif_databases/"
	echo "";
	echo "[TF_list] = A file with a single line of comma-separated Transcription Factor gene symbols. NOTE: truncate with ','";
	echo "";
	exit 1;
fi

dVCF="/research/projects/nicholsgrp/FamilialStudyMethod/common/SOP/JRM_files/fromOriginalVCF_alias_062220/sjhl.sed.hg38_multianno.slivar05.vcf.gz"
dREFERENCE="/research/rgs01/reference/public/genomes/Homo_sapiens/GRCh38/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
dMOTIFDB="/hpcf/authorized_apps/rhel7_apps/meme/install/5.1.0/db/motif_databases/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme"
dTFFILE="/research/projects/cab/automapper/common/jmyers3/FIMOsnp/newVersion/TFs.txt"
dOUTPUT="output_FIMOsnp.txt"

INPUT=$1
OUTPUT=${2:-$dOUTPUT}
VCF=${3:-$dVCF}
REFERENCE=${4:-$dREFERENCE}
MOTIFDB=${5:-$dMOTIFDB}
TFFILE=${6:-$dTFFILE}

TFLIST=$(<${TFFILE})

# LENGTH OF PADDING +/- AROUND VARIANT
PAD=30

# COLUMN OF THE ENCODE TFBS ANNOTATION -1
ENCTFBSPOS=182

# LOAD REQUIRE MODULES
module purge
module load tabix
module load samtools/1.10
module load bcftools/1.10.2
module load meme/5.1.0

sed 1d ${INPUT} | while read -r line; do
	read -ra CURVAR <<< "$line"

	# VARIANT INPUT
	CHROM=${CURVAR[0]}
	POS=${CURVAR[1]}
	REF=${CURVAR[4]}
	ALT=${CURVAR[5]}
	TFBS=${CURVAR[$ENCTFBSPOS]}

	if [ ${#REF} -gt ${#ALT} ]
	then
		OFFSET=`expr ${#REF} + ${PAD}`
	else
		OFFSET=`expr ${#ALT} + ${PAD} - 1`
	fi
	
	# SET REGION TO SCREEN
	START=`expr ${POS} - ${PAD}`
	END=`expr ${POS} + ${OFFSET}`
	REGION="${CHROM}:${START}-${END}"
	echo "Testing Region: ${REGION}"
	
	# SET VARIANT LOCATION
	VAR="${CHROM}:${POS}-${POS}"
	
	# SET FILE NAMES
	VARVCF="${CHROM}_${POS}-${POS}.vcf.gz"
	FORFIMO="${CHROM}_${START}-${END}.fasta"
	FIMORESULTS="fimosnp_${CHROM}-${POS}"
	
	# ISOLATE VARIANT FROM VCF
	tabix -h -p vcf ${VCF} ${VAR} | bgzip >${VARVCF}
	tabix -p vcf ${VARVCF}
	
	# CREATE FASTA FOR REGION WITH REF AND ALT ALLELES
	samtools faidx ${REFERENCE} ${REGION} > ${FORFIMO}
	samtools faidx ${REFERENCE} ${REGION} | bcftools consensus -p var_ ${VARVCF} >> ${FORFIMO} 2> /dev/null
	
	# RUN FIMO
	fimo --verbosity 1 -o ${FIMORESULTS} ${MOTIFDB} ${FORFIMO}

	sort -nk9,9 ${FIMORESULTS}/fimo.tsv > ${FIMORESULTS}/srt_fimo.tsv
	
	# GET HANDLE ON TOP HIT FOR EACH SEQUENCE
	TOPREF=`grep -P "\t${REGION}\t" ${FIMORESULTS}/srt_fimo.tsv -m1`
	TOPALT=`grep -P "\tvar_${REGION}\t" ${FIMORESULTS}/srt_fimo.tsv -m1`
	read -ra TOPREFsplit <<< "${TOPREF}"
	read -ra TOPALTsplit <<< "${TOPALT}"
	
	UPTOPREF=${TOPREFsplit[1]^^}
	UPTOPALT=${TOPALTsplit[1]^^}
	
	# COMPARE THE TOP HIT FOR REF AND ALT THEN OUTPUT YES/NO IF THEY ARE DIFFERENT
	if [ "${TOPREFsplit[0]}" != "${TOPALTsplit[0]}" ]
	then
		if [ -z "${TOPALTsplit[8]}" ] || (( $(echo "${TOPREFsplit[8]} < ${TOPALTsplit[8]}" | bc -l 2>/dev/null) ))
		then
			# LOSS
			if [[ "${TFBS}" == *"${UPTOPREF/::*/}_"* ]]
			then
				echo -e "${line}\tLOSS-${UPTOPREF}" >> ${OUTPUT}
			else
				echo -e "${line}\tLOSS" >> ${OUTPUT}
			fi
		else
			# GAIN
			if [[ "${TFLIST}" == *"${UPTOPALT/::*/},"* ]]
			then
				echo -e "${line}\tGAIN-${UPTOPALT}" >> ${OUTPUT}
			else
				echo -e "${line}\tGAIN" >> ${OUTPUT}
			fi
		fi
	else
		REFCOUNT=`grep -P "\t${REGION}\t" ${FIMORESULTS}/srt_fimo.tsv -c`
		ALTCOUNT=`grep -P "\tvar_${REGION}\t" ${FIMORESULTS}/srt_fimo.tsv -c`

		if [ ${REFCOUNT} -eq ${ALTCOUNT} ]
		then
			echo -e "${line}\tNOCHANGE" >> ${OUTPUT}
		else
			echo -e "${line}\tCOMPLEX" >> ${OUTPUT}
		fi
	fi

	mv ${FORFIMO} ${VARVCF}* ${FIMORESULTS}
	
done 
	
exit 0
