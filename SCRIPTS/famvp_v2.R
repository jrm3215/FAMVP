library(dplyr)
library(readr)
library(kinship2)

args = commandArgs(trailingOnly=TRUE)

### HANDLE ARGUMENTS
if(length(args)!=2){
	stop("Must provide a familyID and a cellType to proceed!")
} else {
	famID <- args[1]
	cellType <- args[2]
}

### READ DATA
rawData <-read_tsv(paste0(famID,".slivar.hg38_multianno.tab.txt"), col_names = TRUE, trim_ws = TRUE, skip_empty_rows = TRUE, na=c("",".","NA","na"), guess_max=100000)
rawData <- rawData %>% filter(is.na(OLD_MULTIALLELIC))

### READ DATABASES AND GENE LISTS
dbLocation <- "/rgs01/project_space/nicholsgrp/FamilialStudyMethod/common/SOP/DATABASES/"
HLfile<- paste0(dbLocation, "HL_candidateGenes_072220.txt")
HLgenes <- read_tsv(HLfile, col_names = TRUE, trim_ws = TRUE, skip_empty_rows = TRUE, na=c("",".","NA","na"))
SJfile <- paste0(dbLocation, "SJFAMILY_Gene_List_160.txt")
SJgenes <- read_tsv(SJfile, col_names = TRUE, trim_ws = TRUE, skip_empty_rows = TRUE, na=c("",".","NA","na"))
gnConst <- read_tsv("/rgs01/project_space/cab/Control/common/reference/gnomAD/gnomad_constraint.tab", col_names = TRUE, trim_ws = TRUE, skip_empty_rows = TRUE, na=c("",".","NA","na"))
gnConst <- gnConst %>% rename(Gene.refGene=gene)

### COMBINE DATA AND DATABASES/GENELISTS
combined <- rawData %>% left_join(HLgenes, by='Gene.refGene') %>% left_join(SJgenes, by='Gene.refGene') %>% left_join(gnConst, by='Gene.refGene')
combined <-  distinct(combined, CHROM, POS,ID,REF,ALT, .keep_all=TRUE)
rm(rawData)

### ADD COLUMNS
combined <- mutate(combined, REF_length=0)
combined <- mutate(combined, ALT_length=0)
combined <- mutate(combined, var_type="none")
combined <- mutate(combined, highRegBase=ifelse(REG_PHRED>10 | CAN_PHRED>10 | PAT_PHRED>10,"yes",NA))
combined <- mutate(combined, fixedGene=NA)
combined <- mutate(combined, encTFBScand=NA)

combined <- mutate(combined, affAC=0)
combined <- mutate(combined, unaffAC=0)

combined <- mutate(combined, famAC=0)
combined <- mutate(combined, famCarrierAC=0)
combined <- mutate(combined, famAF=0)

combined <- mutate(combined, sibAC=0)
combined <- mutate(combined, sibCarrierAC=0)
combined <- mutate(combined, sibAF=0)
combined <- mutate(combined, sibEP=0)
combined <- mutate(combined, sibRM="no")

combined <- mutate(combined, otherAC=0)
combined <- mutate(combined, otherCarrierAC=0)
combined <- mutate(combined, otherAF=0)

combined <- mutate(combined, Priority=5)
combined <- mutate(combined, Rationale="lowPriority")
combined <- mutate(combined, Priority_nc=9)
combined <- mutate(combined, Rationale_nc="lowPriority")


### USEFUL LISTS
LOF <- c("frameshift_insertion", "frameshift_deletion","stopgain")
missense <- c("nonsynonymous_SNV")
INDEL <- c("nonframeshift_insertion","nonframeshift_deletion")
nonDel <- c("synonymous_SNV")
exonicList <-c("exonic","exonic\\x3bsplicing","splicing")
ncrnaList <- c("ncRNA_exonic","ncRNA_exonic\\x3bsplicing","ncRNA_splicing")
ncList <- c("UTR3", "UTR5","downstream","intergenic","intronic","ncRNA_exonic","ncRNA_intronic","ncRNA_splicing","upstream")


############ CODING PRIORITIZATION #############
codingPriority <- function(i){

	combined$fixedGene[i] <- cur[1]

	if(combined$fixedGene[i] %in% HLgenes$Gene.refGene){
		combined$HL_candidate[i] <- "yes"
	}
	if(combined$fixedGene[i] %in% SJgenes$Gene.refGene){
		combined$list160[i] <- "yes"
	}
	if(is.na(combined$REVEL_score[i])){
	   combined$REVEL_score[i] <- 0
	}
	if(is.na(combined$CADD_phred[i])){
	   combined$CADD_phred[i] <- 0
	}

	if((combined$ExonicFunc.refGene[i] %in% c("nonsynonymous_SNV","frameshift_insertion", "frameshift_deletion","stopgain") | combined$Func.refGene[i]=='splicing') &
	   (combined$REVEL_score[i]>=0.5 |  combined$CADD_phred[i]>=20) &
	   (!is.na(combined$HL_candidate[i]) | !is.na(combined$list160[i]))){
			combined$Priority[i] <- 1
			combined$Rationale[i] <- "HL.SJFAM.LOF.INDEL.Nonsyn.REVEL.CADD"
	} else if((combined$ExonicFunc.refGene[i] %in% c("nonsynonymous_SNV","frameshift_insertion", "frameshift_deletion","stopgain") | combined$Func.refGene[i]=='splicing') &
			  (!is.na(combined$HL_candidate[i]) | !is.na(combined$list160[i]))){
			combined$Priority[i] <- 2
			combined$Rationale[i] <- "HL.SJFAM.LOF.INDEL.Nonsyn"
	} else if(combined$ExonicFunc.refGene[i] %in% c("frameshift_insertion", "frameshift_deletion","stopgain") | combined$Func.refGene[i]=='splicing'){
			combined$Priority[i] <- 3
			combined$Rationale[i] <- "LOF.INDEL"
	} else if(combined$ExonicFunc.refGene[i] %in% c("nonsynonymous_SNV") & (combined$REVEL_score[i]>=0.5 | combined$CADD_phred[i]>=20)){
			combined$Priority[i] <- 4
			combined$Rationale[i] <- "Nonsyn.REVEL.CADD"
	}

	return(combined)
}


################ NON_CODING PRIORITIZATION ##############
nonCodingPriority <- function(i){

	if(combined$Func.refGene[i]=="intergenic"){
		checkDist <- unlist(strsplit(gsub("dist\\\\x3d","",combined$GeneDetail.refGene[i]), "\\\\x3b"))
		if(cur[1]=='NONE' & cur[2]=='NONE'){
			combined$fixedGene[i] <- cur[1]
		}else if(cur[1]=='NONE'){
			combined$fixedGene[i] <- cur[2]
		}else if(cur[2]=='NONE'){
			combined$fixedGene[i] <- cur[1]
		}else if(as.numeric(checkDist[1]) < as.numeric(checkDist[2])){
			combined$fixedGene[i] <- cur[1]
		} else{
			combined$fixedGene[i] <- cur[2]
		}
	} else {
		combined$fixedGene[i] <- cur[1]
	}

	if(combined$fixedGene[i] %in% HLgenes$Gene.refGene){
		combined$HL_candidate[i] <- "yes"
	}
	if(combined$fixedGene[i] %in% SJgenes$Gene.refGene){
		combined$list160[i] <- "yes"
	}

	if((!is.na(combined$wgEncodeRegDnaseClustered[i]) | !is.na(combined$encRegTfbsClustered[i])) & (!is.na(combined$highRegBase[i]) | combined$var_type[i]=='INDEL') & (grepl(cellType, combined$wgEncodeRegDnaseClustered[i]) | grepl(cellType, combined$encRegTfbsClustered[i])) & (!is.na(combined$HL_candidate[i]) | !is.na(combined$list160[i]))){
		combined$Priority_nc[i] <- 1
		combined$Rationale_nc[i] <- "Regulatory_highImpact_CellType_KnownGene"
	} else if((!is.na(combined$wgEncodeRegDnaseClustered[i]) | !is.na(combined$encRegTfbsClustered[i])) & (!is.na(combined$highRegBase[i]) | combined$var_type[i]=='INDEL') & (!is.na(combined$HL_candidate[i]) | !is.na(combined$list160[i]))){
		combined$Priority_nc[i] <- 2
		combined$Rationale_nc[i] <- "Regulatory_highImpact_KnownGene"
	} else if((!is.na(combined$wgEncodeRegDnaseClustered[i]) | !is.na(combined$encRegTfbsClustered[i])) & (!is.na(combined$highRegBase[i]) | combined$var_type[i]=='INDEL') & (grepl(cellType, combined$wgEncodeRegDnaseClustered[i]) | grepl(cellType, combined$encRegTfbsClustered[i]))){
		combined$Priority_nc[i] <- 3
		combined$Rationale_nc[i] <- "Regulatory_highImpact_CellType"
	} else if((!is.na(combined$wgEncodeRegDnaseClustered[i]) | !is.na(combined$encRegTfbsClustered[i])) & (!is.na(combined$HL_candidate[i]) | !is.na(combined$list160[i]))){
		combined$Priority_nc[i] <- 4
		combined$Rationale_nc[i] <- "Regulatory_KnownGene"
	} else if((!is.na(combined$wgEncodeRegDnaseClustered[i]) | !is.na(combined$encRegTfbsClustered[i])) & (!is.na(combined$highRegBase[i]) | combined$var_type[i]=='INDEL')){
		combined$Priority_nc[i] <- 5
		combined$Rationale_nc[i] <- "Regulatory_highImpact"
	} else if((!is.na(combined$wgEncodeRegDnaseClustered[i]) | !is.na(combined$encRegTfbsClustered[i])) & (grepl(cellType, combined$wgEncodeRegDnaseClustered[i]) | grepl(cellType, combined$encRegTfbsClustered[i]))){
		combined$Priority_nc[i] <- 6
		combined$Rationale_nc[i] <- "Regulatory_CellType"
	} else if((!is.na(combined$wgEncodeRegDnaseClustered[i]) | !is.na(combined$encRegTfbsClustered[i]))){
		combined$Priority_nc[i] <- 7
		combined$Rationale_nc[i] <- "Regulatory"
	} else if(combined$Func.refGene[i] %in% ncrnaList){
		combined$Priority_nc[i] <- 8
		combined$Rationale_nc[i] <- "ncRNA_overlap"
	} 

	if(combined$fixedGene[i]=='NONE'){
		combined$Priority_nc[i] <- 9
		combined$Rationale_nc[i] <- "low"
	}

	if(!is.na(combined$encRegTfbsClustered[i])){
		for(t in 1:nrow(HLgenes)){
			if(grepl(paste0(HLgenes$Gene.refGene[t],"_"),combined$encRegTfbsClustered[i])){
				if(!is.na(combined$encTFBScand[i])){
					combined$encTFBScand[i] <- paste0(combined$encTFBScand[i],",",HLgenes$Gene.refGene[t])
				} else {
					combined$encTFBScand[i] <- HLgenes$Gene.refGene[t]
				}
			}
		}
	}

	return(combined)
}



######## PENETRANCE ############

pedFile	<- paste0(famID, ".ped")
aliasFile	<- paste0(famID, ".alias")

aliasData = read.delim(aliasFile, header = T, as.is = T, sep = "\t")
if(is.na(aliasData$X.affecteds)){ aliasData$X.affecteds <- ''}
if(is.na(aliasData$unaffecteds)){ aliasData$unaffecteds <- ''}
affecteds <-unlist(strsplit(aliasData$X.affecteds, ","))
unaffecteds <-unlist(strsplit(aliasData$unaffecteds, ","))

pedData = read.table(pedFile, header = T, as.is = T, sep = "\t")
test <- pedigree(id=pedData$Individual.ID, dadid=pedData$PaternalID, momid=pedData$MaternalID, sex=pedData$Gender, famid=pedData$Family.ID, affected=pedData$Phenotype,missid=0)
IDs <-test[famID]$id
IDs

affectedData <- subset(pedData, Phenotype==2)
unaffectedData <- subset(pedData, Phenotype==1)
siblings <- vector()
for(t in 1:nrow(unaffectedData)){
    if(unaffectedData$PaternalID[t] %in% affectedData$PaternalID | unaffectedData$MaternalID[t] %in% affectedData$MaternalID){
        siblings <- c(siblings, unaffectedData$Individual.ID[t])
    }
}

### START LOOP ###
famAC <- function(i) {

	case <- slice(combined,i)
	findIDs <-gsub("$",".GT",IDs)
	gtCase <- select(case, one_of(findIDs))
	foundGTs <- unlist(gtCase)
	alleles <- select(case, REF, ALT) %>% unlist(use.names=FALSE)

	if(alleles[2]=="<DEL>"){
		alleles[2]<-""
	}

	het <- paste(alleles, collapse="/")
	alt <- paste(alleles[2], alleles[2], sep="/")
	ref <- paste(alleles[1], alleles[1], sep="/")
	df.gtCase <- as.data.frame(foundGTs,stringsAsFactors=FALSE)
	gtRN <- rownames(df.gtCase)
	rownames(df.gtCase) <-gsub(".GT","",gtRN)
	n = length(df.gtCase$foundGTs)

	for(j in 1:n){ 
		if(df.gtCase$foundGTs[j] == het | df.gtCase$foundGTs[j] == alt){
			df.gtCase$carrier[j] <- 2
		} else if(df.gtCase$foundGTs[j] == ref){
			df.gtCase$carrier[j] <- 1
		} else {
			df.gtCase$carrier[j] <- 0
		}

		if(rownames(df.gtCase[j,]) %in% affecteds){
			df.gtCase$category[j] <- "aff"
		} else if(rownames(df.gtCase[j,]) %in% unaffecteds){
			df.gtCase$category[j] <- "unaff"
		} else if(rownames(df.gtCase[j,]) %in% siblings){
			df.gtCase$category[j] <- "sib"
		} else {
			df.gtCase$category[j] <- "other"
		}
	}

	#### counts
	combined$famAC[i] <- length(df.gtCase[which(df.gtCase$carrier>0),]$category)
	combined$affAC[i] <- length(df.gtCase[which(df.gtCase$category=='aff' & df.gtCase$carrier>0),]$category)
	combined$unaffAC[i] <- length(df.gtCase[which(df.gtCase$category=='unaff' & df.gtCase$carrier>0),]$category)
	combined$sibAC[i] <- length(df.gtCase[which(df.gtCase$category=='sib' & df.gtCase$carrier>0),]$category)
	combined$otherAC[i] <- length(df.gtCase[which(df.gtCase$category=='other' & df.gtCase$carrier>0),]$category)
	combined$famCarrierAC[i] <- length(df.gtCase[which(df.gtCase$carrier==2),]$carrier)
	combined$sibCarrierAC[i] <- length(df.gtCase[which(df.gtCase$category=='sib' & df.gtCase$carrier==2),]$carrier)
	combined$otherCarrierAC[i] <- length(df.gtCase[which(df.gtCase$category=='other' & df.gtCase$carrier==2),]$carrier)

	#### family allele frequency
	combined$famAF[i] <- combined$famCarrierAC[i] / combined$famAC[i]

	#### sibling allele frequency
	if(combined$sibAC[i]>0){
		combined$sibAF[i] <- combined$sibCarrierAC[i] / combined$sibAC[i]
	} else {
		combined$sibAF[i] <- 0
	}

	#### other family allele frequency
	if(combined$otherAC[i]>0){
		combined$otherAF[i] <- combined$otherCarrierAC[i] / combined$otherAC[i]
	} else {
		combined$otherAF[i] <- 0
	}

	#### Estimated penetrance based on sibiling carrier frequency
	if(combined$sibAF[i]==1){
		combined$sibEP[i] <- (-100)
	} else if(combined$sibAF[i]==0.5){
		combined$sibEP[i] <- 0
	} else {
		combined$sibEP[i] <- (((2 * combined$sibAF[i]) - 1) / (combined$sibAF[i] - 1))
	}

	if(combined$sibAC[i] > 1 & combined$sibEP[i] < 0){
		combined$sibRM[i] <- "yes"
	}

	return(combined)
}

### MAIN LOOP
for(i in 1:nrow(combined)){
	print(i)
	combined$REF_length[i] <- nchar(combined$REF[i])
	combined$ALT_length[i] <- nchar(combined$ALT[i])
	if(combined$REF_length[i] ==1 & combined$ALT_length[i] ==1){
		combined$var_type[i]<-"SNP"
	} else {
		combined$var_type[i]<-"INDEL"
	}

	cur <- unlist(strsplit(combined$Gene.refGene[i], "\\\\x3b"))

	if(combined$Func.refGene[i] %in% exonicList){
		(combined <- codingPriority(i))
	} else {
		(combined <- nonCodingPriority(i))
	}

	(combined <- famAC(i))
}

### WRITE PRIORITIZED VARIANT FILES
priority.coding <- combined %>% filter(Func.refGene %in% exonicList & Priority < 5) %>% filter(is.na(gnomad211_exome_AF) | (gnomad211_exome_AF < 0.01 & gnomad211_exome_AF_afr < 0.01 & gnomad211_exome_AF_nfe < 0.01 & gnomad211_exome_AF_amr < 0.05 & gnomad211_exome_AF_asj < 0.05 & gnomad211_exome_AF_eas < 0.05 & gnomad211_exome_AF_fin < 0.05 & gnomad211_exome_AF_oth < 0.05 & gnomad211_exome_AF_sas < 0.05))
codingStr <- paste0(famID,"_priority_coding.txt")
write_tsv(priority.coding, codingStr)

print("Coding Priority Table")
print(table(priority.coding$Priority))

print("Coding variants to remove because of estimated penetrance:")
print(nrow(priority.coding %>% filter(sibRM=="yes")))

if(nrow(priority.coding) >1){
	pdf(paste0(famID,"_priority_coding_EP.pdf"))
	cEP <- density(priority.coding$sibEP)
	plot(cEP, main=paste0(famID,": Priority Coding Variant Estimated Penetrance"))
	dev.off()
}

priority.noncoding <- combined %>% filter(!(Func.refGene %in% exonicList) & Priority_nc < 9) %>% filter(is.na(gnomad30_genome_AF) | (gnomad30_genome_AF < 0.01 & gnomad30_genome_AF_afr < 0.01 & gnomad30_genome_AF_nfe < 0.01 & gnomad30_genome_AF_ami < 0.05 & gnomad30_genome_AF_amr < 0.05 & gnomad30_genome_AF_asj < 0.05 & gnomad30_genome_AF_eas < 0.05 & gnomad30_genome_AF_fin < 0.05 & gnomad30_genome_AF_oth < 0.05 & gnomad30_genome_AF_sas < 0.05))
noncodingStr <- paste0(famID,"_priority_noncoding.txt")
write_tsv(priority.noncoding, noncodingStr)

print("Non-Coding Priority Table")
print(table(priority.noncoding$Priority_nc))

print("Non-Coding variants to remove because of estimated penetrance:")
print(nrow(priority.noncoding %>% filter(sibRM=="yes")))

if(nrow(priority.noncoding) >1){
	pdf(paste0(famID,"_priority_noncoding_EP.pdf"))
	ncEP <- density(priority.noncoding$sibEP)
	plot(ncEP, main=paste0(famID,": Priority Non-Coding Variant Estimated Penetrance"))
	dev.off()
}

### WRITE R DATA OBJECT
imageStr <- paste0(famID, "_prioritized.Rdata")
save.image(imageStr)
