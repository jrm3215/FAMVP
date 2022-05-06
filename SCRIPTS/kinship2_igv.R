library(dplyr)
library(readr)
library(stringr)
library(kinship2)

args = commandArgs(trailingOnly=TRUE)

fam <- basename(getwd())
priorityFile <- args[1]
pedFile	<- paste0(fam, ".ped")
aliasFile	<- paste0(fam, ".alias")

aliasData = read.delim(aliasFile, header = T, as.is = T, sep = "\t")
if(is.na(aliasData$X.affecteds)){ aliasData$X.affecteds <- ''}
if(is.na(aliasData$unaffecteds)){ aliasData$unaffecteds <- ''}
affecteds <-unlist(strsplit(aliasData$X.affecteds, ","))
unaffecteds <-unlist(strsplit(aliasData$unaffecteds, ","))

pedData = read.table(pedFile, header = T, as.is = T, sep = "\t")
test <- pedigree(id=pedData$Individual.ID, dadid=pedData$PaternalID, momid=pedData$MaternalID, sex=pedData$Gender, famid=pedData$Family.ID, affected=pedData$Phenotype,missid=0)
IDs <-test[fam]$id
IDs
IDs1 <-gsub("_G.-.*","",IDs)
IDs1

### READ IN THE PRIORITIZED VARIANTS ###
tabFile <- read_tsv(priorityFile, col_names = TRUE, trim_ws = TRUE, skip_empty_rows = TRUE, na=c("",".","NA","na"), guess_max=100000)
if(args[2] == "coding"){
	outputDIR <-paste0(fam,"_coding_kinship2")
	dir.create(outputDIR)
	setwd(outputDIR)
	subTabFile <- tabFile %>% filter(Priority <= 3)
	if(nrow(subTabFile) == 0) {stop("There are no coding variants with priority < 3")}
	rm(tabFile)
	highPriority <- gsub("priority_coding.txt","highPriority_coding.txt",priorityFile)
	write_tsv(subTabFile, highPriority)
} else if(args[2] == "noncoding"){
	outputDIR <-paste0(fam,"_noncoding_kinship2")
	dir.create(outputDIR)
	setwd(outputDIR)
	subTabFile <- tabFile %>% filter(Priority_nc < 3)
	if(nrow(subTabFile) == 0) {stop("There are no non-coding variants with priority < 3")}
	rm(tabFile)
	highPriority <- gsub("priority_noncoding.txt","highPriority_noncoding.txt",priorityFile)
	write_tsv(subTabFile, highPriority)
	fimoRes <- gsub(".txt","_fimo.txt",highPriority)
	fimoCMD <- paste0("../../SCRIPTS/FIMOsnp.sh ", highPriority, " ", fimoRes)
	system(fimoCMD)
}

### IF A BAM LIST IS GIVEN CREATE IGV SNAPSHOTS ###
if(length(args) == 3){
	bams<- read_tsv(args[3], col_names=F)
	write("new",file="IGV_snapshots.bat")
	write("genome hg38",file="IGV_snapshots.bat",append=TRUE)
	write("maxPanelHeight 100",file="IGV_snapshots.bat",append=TRUE)
	write("snapshotDirectory .",file="IGV_snapshots.bat",append=TRUE)

	for(curSamp in IDs1){
		curBAM<-filter(bams, str_detect(X1, curSamp)) %>% select(X1) %>% unlist(use.names=F)
		write(paste0("load ",curBAM),file="IGV_snapshots.bat",append=TRUE)
	}
}

findIDs <-gsub("$",".GT",IDs)
### START LOOP ###
for(iter in 1:nrow(subTabFile)) {
	case <- slice(subTabFile,iter)
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
	n = length(df.gtCase$foundGTs)

	for(i in 1:n){ 
		if(df.gtCase$foundGTs[i] == ref){ 
			   df.gtCase$hetCode[i] <- 0
			   df.gtCase$altCode[i] <- 0
			   df.gtCase$misCode[i] <- 0
		} else if(df.gtCase$foundGTs[i] == het){
			   df.gtCase$hetCode[i] <- 1
			   df.gtCase$altCode[i] <- 0
			   df.gtCase$misCode[i] <- 0
		} else if(df.gtCase$foundGTs[i] == alt){
			   df.gtCase$hetCode[i] <- 0
			   df.gtCase$altCode[i] <- 1
			   df.gtCase$misCode[i] <- 0
		} else {
			   df.gtCase$hetCode[i] <- 0
			   df.gtCase$altCode[i] <- 0
			   df.gtCase$misCode[i] <- 1
		}
	}

	df.aff <- data.frame(affected=test[fam]$affected)
	rownames(df.aff) <- test[fam]$id
	gtRN <- rownames(df.gtCase)
	rownames(df.gtCase) <-gsub(".GT","",gtRN)
	merged <- merge(df.aff, df.gtCase, by=0,all=TRUE,sort=FALSE)
	srt_merged <- merged[match(rownames(df.aff), merged$Row.names),]
	srt_merged[is.na(srt_merged)] <- 0
	sub_srt_merged <- srt_merged[c("affected","hetCode","altCode","misCode")]

	j = length(srt_merged$Row.names)
	for(k in 1:j){ 
		if(srt_merged$Row.names[k] %in% affecteds){ 
			   srt_merged$sampColors[k] <- "#63ACBE"
		} else if(srt_merged$Row.names[k] %in% unaffecteds){
			   srt_merged$sampColors[k] <- "#EE442F"
		} else {
			   srt_merged$sampColors[k] <- "#000000"
		}
	}

	sampColors <- srt_merged$sampColors

	### CREATE CASE PEDIGREE ###
	current <- pedigree(id=pedData$Individual.ID, dadid=pedData$PaternalID, momid=pedData$MaternalID, sex=pedData$Gender, famid=pedData$Family.ID, affected=as.matrix(sub_srt_merged),missid=0)

	### PLOT PEDIGREE ###
	caseAttr <- select(case, CHROM, POS, REF, ALT) %>% unlist(use.names=FALSE)
	plotFile <- paste(fam, caseAttr[1], caseAttr[2], caseAttr[3], caseAttr[4],  "PED.png", sep="_")

	png(plotFile, width = 11, height = 8, units = 'in',res = 600, type = "cairo")
	plot.pedigree(current[fam],cex = 1, id=IDs1,mar = c(2.1, 2.1, 2.1, 2.1),symbolsize=3, col=sampColors)
	pedigree.legend(current, location="topright", radius=.1,labels=c("Affected","HET","ALT","MISSING"))
	dev.off()

	if(length(args) == 3){
		write(paste0("goto ", caseAttr[1], ":", caseAttr[2], "-", caseAttr[2]),file="IGV_snapshots.bat",append=TRUE)
		write("collapse",file="IGV_snapshots.bat",append=TRUE)
		snapFile <- gsub("PED.png","IGV.png",plotFile)
		write(paste0("snapshot ", snapFile),file="IGV_snapshots.bat",append=TRUE)
	}

}

### CREATE IGV SNAPSHOTS
if(length(args) == 3){
	write("exit",file="IGV_snapshots.bat",append=TRUE)
	system('xvfb-run -s "-screen 0 1280x720x8" -d java -Xmx4000m -jar /hpcf/authorized_apps/rhel7_apps/igv/install/2.3.2/igv.jar -b IGV_snapshots.bat')
}


### END ###
