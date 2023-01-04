suppressMessages(library(plyr))
suppressMessages(library(PharmacoGx))
#library(matrixStats)
options(stringsAsFactors = FALSE)


###########################################################
##### Determine biological replicates across datasets #####
###########################################################

#match cell.id to unique.cellid from PharmacoGx
matchToIDTable <- function(ids,tbl, column, returnColumn="unique.cellid") {
  sapply(ids, function(x) {
    myx <- grep(paste0("((///)|^)",Hmisc::escapeRegex(x),"((///)|$)"), tbl[,column])
    if(length(myx) > 1){
      stop("Something went wrong in curating ids, we have multiple matches")
    }
    if(length(myx) == 0){return(NA_character_)}
    return(tbl[myx, returnColumn])
  })
}


#read in current cell annotations
cell_all <- read.csv(file = "/data/rnaseq_meta/cell_annotation_all.csv", na.strings=c("", " ", "NA"))

#read in gcsi cell annotations
gcsi <- read.csv(file = "/data/rnaseq_meta/gcsi_rnaseq_meta.csv")
gcsi$cellid <- matchToIDTable(ids=gcsi$Cell_line , tbl = cell_all, column = "GNE.cellid", returnColumn = "unique.cellid")
rownames(gcsi) <- gcsi$alias

#read in ccle cell annotations
ccle <- read.csv(file = "/data/rnaseq_meta/ccle_rnaseq_meta.csv")
ccle$cellid <- matchToIDTable(ids=ccle$Cell_Line , tbl = cell_all, column = "CCLE.cellid", returnColumn = "unique.cellid")
rownames(ccle) <- ccle$Run

#read in gdsc cell annotations
gdsc <- read.csv(file = "/data/rnaseq_meta/gdsc_rnaseq_meta.txt", sep = "\t")
gdsc <- gdsc[which(!gdsc$Comment.SUBMITTED_FILE_NAME. == "15552_5.cram"),]
gdsc$cellid <- matchToIDTable(ids=gdsc$Source.Name, tbl = cell_all, column = "GDSC_rnaseq.cellid", returnColumn = "unique.cellid")
gdsc$files <- gsub(".cram","",gdsc$Comment.SUBMITTED_FILE_NAME.)
rownames(gdsc) <- gdsc$files

#filter metadata -> keep intersected cell lines only
intersected_rnacells <- Reduce(intersect, list(gdsc$cellid, ccle$cellid, gcsi$cellid))

gcsi <- gcsi[which(gcsi$cellid %in% intersected_rnacells),]
ccle <- ccle[which(ccle$cellid %in% intersected_rnacells),]
gdsc <- gdsc[which(gdsc$cellid %in% intersected_rnacells),]

save(intersected_rnacells, file="../results/inter.RData") #save biological replicates


########################################################
########## Get circRNA counts from CIRI2 data ##########
########################################################

circRNA_counts_CIRI <- function(dir_path){
  
  sample_names <- c() #holds all sample names
  circRNA_counts_df <- data.frame(matrix(ncol=0, nrow = 0))
  
  ciri_files <- list.files(path=dir_path,  pattern = "\\.tsv$", full.names = T) 
  for (f in 1:length(ciri_files)) {
    sample <- read.table(file = ciri_files[f], 
                         sep = '\t', 
                         skip = 1,
                         header = FALSE)
    # store sample names
    #sample_name <- gsub("\\..*","", ciri_files[f])
    sample_name <- gsub(".*/","", ciri_files[f])
    sample_name <- gsub("*.tsv", "", sample_name)
    sample_names <- c(sample_names, sample_name)
    
    # keep only exonic circRNAs and circRNAs with >=2 junction reads
    colnames(sample) <- c("circRNA_ID","chr","circRNA_start","circRNA_end","junction_reads", "SM_MS_SMS", "non_junction_reads", "junction_reads_ratio", "circRNA_type", "gene_id", "strand","junction_reads_ID")
    sample <- sample[which(sample$circRNA_type=="exon"),] 
    sample <- sample[which(sample$junction_reads >= 2),] 
    
    #convert coordinate from 1-based to 0-based to match CIRCexplorer2 coordinates
    sample$circRNA_start <- sample$circRNA_start - 1
    sample$circRNA_ID <- paste0(sample$chr, ":", sample$circRNA_start, "|", sample$circRNA_end)
    
    #keep only # junction reads (counts) and circRNA IDs
    sample <- sample[,c("circRNA_ID", "junction_reads")]
    sample <- t(sample)
    colnames(sample) <- sample[1,]
    sample <- as.data.frame(sample)
    sample <- sample[-c(1),]
    
    circRNA_counts_df <- rbind.fill(circRNA_counts_df, sample)
    
    
  }
  circRNA_counts_df <- as.data.frame(circRNA_counts_df)
  circRNA_counts_df[] <- sapply(circRNA_counts_df, as.integer)
  
  
  rownames(circRNA_counts_df) <- sample_names 
  return(circRNA_counts_df)
} 

ciri_gcsi <- circRNA_counts_CIRI(dir_path = "/data/processed_cellline/CIRI2/gCSI/result")
ciri_ccle <- circRNA_counts_CIRI(dir_path = "/data/processed_cellline/CIRI2/CCLE/result")
ciri_gdsc <- circRNA_counts_CIRI(dir_path = "/data/processed_cellline/CIRI2/GDSC/result")

ciri_gcsi[is.na(ciri_gcsi)] <- 0
ciri_ccle[is.na(ciri_ccle)] <- 0
ciri_gdsc[is.na(ciri_gdsc)] <- 0

#clean up sample names (remove 'gcsi' & 'gdsc' strings from samples)
rownames(ciri_gcsi) <- gsub("gcsi", "", rownames(ciri_gcsi))
rownames(ciri_gdsc) <- gsub("gdsc", "", rownames(ciri_gdsc))


#average SR technical replicate circRNA count in gCSI
ciri_gcsi$cellid <- gcsi$cellid[match(rownames(ciri_gcsi), rownames(gcsi))]

gcsi_SR_mean <- c() #hold means for each circRNA transcript
for (k in 1:(length(colnames(ciri_gcsi)) -1)) {
  mean_count <- mean(as.numeric(ciri_gcsi[which(ciri_gcsi$cellid=="SR"),k]))
  gcsi_SR_mean <- c(gcsi_SR_mean, mean_count)
}
gcsi_SR_mean <- c(gcsi_SR_mean, "SR")

# replace two SR rows with mean
ciri_gcsi <- ciri_gcsi[-which(ciri_gcsi$cellid=="SR"),]
ciri_gcsi <- rbind(ciri_gcsi, gcsi_SR_mean)
rownames(ciri_gcsi) <- c(rownames(ciri_gcsi)[1:47], "636074")


#get unique.cellids for samples + order CCLE & GDSC
ciri_ccle$cellid <- ccle$cellid[match(rownames(ciri_ccle), rownames(ccle))]
ciri_gdsc$cellid <- gdsc$cellid[match(rownames(ciri_gdsc), rownames(gdsc))]

# make unique cellid = rownames and remove $cellid
rownames(ciri_gcsi) <- ciri_gcsi$cellid
rownames(ciri_ccle) <- ciri_ccle$cellid
rownames(ciri_gdsc) <- ciri_gdsc$cellid

ciri_gcsi <- subset(ciri_gcsi, select = -c(cellid))
ciri_ccle <- subset(ciri_ccle, select = -c(cellid))
ciri_gdsc <- subset(ciri_gdsc, select = -c(cellid))

# transpose dataframes to format for stability calculation
ciri_gcsi <- as.data.frame(t(ciri_gcsi))
ciri_ccle <- as.data.frame(t(ciri_ccle))
ciri_gdsc <- as.data.frame(t(ciri_gdsc))

##only keep transcripts found in all datasets now after filitering for 0 expression above (result: 283 circRNA transcripts)
intersected_t <- Reduce(intersect, list(rownames(ciri_gcsi), rownames(ciri_ccle), rownames(ciri_gdsc)))
ciri_gcsi <- ciri_gcsi[intersected_t,]
ciri_ccle <- ciri_ccle[intersected_t,]
ciri_gdsc <- ciri_gdsc[intersected_t,]

# log2 + 1 normalization
c1 <- sapply(ciri_gcsi, as.numeric)
c2 <- sapply(ciri_ccle, as.numeric)
c3 <- sapply(ciri_gdsc, as.numeric)

rownames(c1) <- rownames(ciri_gcsi)
rownames(c2) <- rownames(ciri_ccle)
rownames(c3) <- rownames(ciri_gdsc)

ciri_gcsi <- log2(c1 + 1)
ciri_ccle <- log2(c2 + 1)
ciri_gdsc <- log2(c3 + 1)


################################################################
########## Get circRNA counts from CIRCexplorer2 data ##########
################################################################


circRNA_counts_CIRCexplorer <- function(dir_path){
  
  sample_names <- c() #holds all sample names
  circRNA_counts_df <- data.frame(matrix(ncol=0, nrow = 0))
  
  circ_files <- list.files(path=dir_path,  pattern = ".txt", full.names = TRUE, recursive = TRUE)
  for (f in 1:length(circ_files)) {
    sample <- read.table(file = circ_files[f], 
                         sep = '\t',
                         header = FALSE)
    # store sample names
    sample_name <- sub(".*/ *(.*?) *.txt.*", "\\1", circ_files[f])
    sample_names <- c(sample_names, sample_name)
    
    # keep only exonic circRNAs and circRNAs with >=2 junction reads
    colnames(sample) <- c("chrom","start","end", "name","score","strand","thickStart","thickEnd","itemRgb","exonCount","exonSizes", "exonOffsets", "readNumber", "circType", "geneName", "isoformName", "index", "flankIntron")
    sample <- sample[which(sample$circType=="circRNA"),] 
    sample <- sample[which(sample$readNumber >= 2),] 
    
    sample$circRNA_ID <- paste0(sample$chrom,":",sample$start,"|",sample$end) #create circRNA_id using chr, start & end coordinates
    
    # keep only #reads (counts) and circRNA IDs
    sample <- sample[,c("readNumber", "circRNA_ID")]
    sample <- t(sample)
    colnames(sample) <- sample[2,]
    sample <- as.data.frame(sample)
    sample <- sample[-c(2),]
    
    
    circRNA_counts_df <- rbind.fill(circRNA_counts_df, sample)
    
  }
  circRNA_counts_df <- as.data.frame(circRNA_counts_df)
  circRNA_counts_df[] <- sapply(circRNA_counts_df, as.integer)
  
  
  rownames(circRNA_counts_df) <- sample_names 
  return(circRNA_counts_df)
} 

circ_gcsi <- circRNA_counts_CIRCexplorer(dir_path = "/data/processed_cellline/CIRCexplorer2/gCSI/annotate")
circ_ccle <- circRNA_counts_CIRCexplorer(dir_path = "/data/processed_cellline/CIRCexplorer2/CCLE/annotate")
circ_gdsc <- circRNA_counts_CIRCexplorer(dir_path = "/data/processed_cellline/CIRCexplorer2/GDSC/annotate")


circ_gcsi[is.na(circ_gcsi)] <- 0
circ_ccle[is.na(circ_ccle)] <- 0
circ_gdsc[is.na(circ_gdsc)] <- 0



#clean up sample names
rownames(circ_gcsi) <- sub(".*gcsi *(.*?) *_.*", "\\1", rownames(circ_gcsi))
rownames(circ_ccle) <- sub("* *(.*?) *_circular.*", "\\1", rownames(circ_ccle))
rownames(circ_gdsc) <- sub(".*gdsc *(.*?) *_circular.*", "\\1", rownames(circ_gdsc))

#average SR technical replicate circRNA count in gCSI
circ_gcsi$cellid <- gcsi$cellid[match(rownames(circ_gcsi), rownames(gcsi))]

gcsi_SR_mean <- c() #hold means for each circRNA transcript
for (k in 1:(length(colnames(circ_gcsi)) -1)) {
  mean_count <- mean(as.numeric(circ_gcsi[which(circ_gcsi$cellid=="SR"),k]))
  gcsi_SR_mean <- c(gcsi_SR_mean, mean_count)
}
gcsi_SR_mean <- c(gcsi_SR_mean, "SR")

# replace two SR rows with mean
circ_gcsi <- circ_gcsi[-which(circ_gcsi$cellid=="SR"),]
circ_gcsi <- rbind(circ_gcsi, gcsi_SR_mean)
rownames(circ_gcsi) <- c(rownames(circ_gcsi)[1:47], "636074")


#get unique.cellids for samples + order CCLE & GDSC
circ_ccle$cellid <- ccle$cellid[match(rownames(circ_ccle), rownames(ccle))]
circ_gdsc$cellid <- gdsc$cellid[match(rownames(circ_gdsc), rownames(gdsc))]


# make unique cellid = rownames and remove $cellid
rownames(circ_gcsi) <- circ_gcsi$cellid
rownames(circ_ccle) <- circ_ccle$cellid
rownames(circ_gdsc) <- circ_gdsc$cellid

circ_gcsi <- subset(circ_gcsi, select = -c(cellid))
circ_ccle <- subset(circ_ccle, select = -c(cellid))
circ_gdsc <- subset(circ_gdsc, select = -c(cellid))


# transpose dataframes to format for stability calculation
circ_gcsi <- as.data.frame(t(circ_gcsi))
circ_ccle <- as.data.frame(t(circ_ccle))
circ_gdsc <- as.data.frame(t(circ_gdsc))

##only keep transcripts found in all datasets now after filitering for 0 expression above (result: 280 circRNA transcripts)
intersected_t <- Reduce(intersect, list(rownames(circ_gcsi), rownames(circ_ccle), rownames(circ_gdsc)))
circ_gcsi <- circ_gcsi[intersected_t,]
circ_ccle <- circ_ccle[intersected_t,]
circ_gdsc <- circ_gdsc[intersected_t,]

# log2 + 1 normalization
c1 <- sapply(circ_gcsi, as.numeric)
c2 <- sapply(circ_ccle, as.numeric)
c3 <- sapply(circ_gdsc, as.numeric)

rownames(c1) <- rownames(circ_gcsi)
rownames(c2) <- rownames(circ_ccle)
rownames(c3) <- rownames(circ_gdsc)

circ_gcsi <- log2(c1 + 1)
circ_ccle <- log2(c2 + 1)
circ_gdsc <- log2(c3 + 1)

save(ciri_gcsi, ciri_ccle, ciri_gdsc, circ_gcsi, circ_ccle, circ_gdsc, 
     file = "../results/circRNAcounts.RData")



########################################
#### Lung circRNA transcript counts ####
########################################

lung_ribozero_ciri <- circRNA_counts_CIRI(dir_path = "/data/processed_lung/CIRI2/RiboZero/result")
lung_polyA_ciri <- circRNA_counts_CIRI(dir_path = "/data/processed_lung/CIRI2/polyA/result")

lung_ribozero_ciri[is.na(lung_ribozero_ciri)] <- 0
lung_polyA_ciri[is.na(lung_polyA_ciri)] <- 0

lung_ribozero_circ <- circRNA_counts_CIRCexplorer(dir_path = "/data/processed_lung/CIRCexplorer2/RiboZero/annotate")
lung_polyA_circ <- circRNA_counts_CIRCexplorer(dir_path = "/data/processed_lung/CIRCexplorer2/polyA/annotate")                          


lung_ribozero_circ[is.na(lung_ribozero_circ)] <- 0
lung_polyA_circ[is.na(lung_polyA_circ)] <- 0                                        


# log2 + 1 normalization
l1 <- sapply(lung_ribozero_ciri, as.numeric)
l2 <- sapply(lung_polyA_ciri, as.numeric)
l3 <- sapply(lung_ribozero_circ, as.numeric)
l4 <- sapply(lung_polyA_circ, as.numeric)

rownames(l1) <- rownames(lung_ribozero_ciri)
rownames(l2) <- rownames(lung_polyA_ciri)
rownames(l3) <- rownames(lung_ribozero_circ)
rownames(l4) <- rownames(lung_polyA_circ)

lung_ribozero_ciri <- log2(l1 + 1)
lung_polyA_ciri <- log2(l2 + 1)
lung_ribozero_circ <- log2(l3 + 1)
lung_polyA_circ <- log2(l4 + 1)


save(lung_ribozero_ciri, lung_polyA_ciri, lung_ribozero_circ, lung_polyA_circ, file = "../results/lungcounts.RData")


#########################################
##### Get isoform transcript counts #####
#########################################

#PSets downloaded from ORCESTRA (orcestra.ca) with the following DOI's and not uploaded to code ocean due to size:
#gCSI DOI: 10.5281/zenodo.3905452
#CCLE DOI: 10.5281/zenodo.3905462
#GDSCv2-8.2 DOI: 10.5281/zenodo.3905481

#read in PSets
gCSI <- readRDS("/data/PSets/gCSI.rds")
CCLE <- readRDS("/data/PSets/CCLE.rds")
GDSC <- readRDS("/data/PSets/GDSC2-8.2.rds")


#summarize gene expression across replicates for each pSet
gCSI_rna <- summarizeMolecularProfiles(gCSI, mDataType = "Kallisto_0.46.1.isoforms.counts", cell.lines = intersected_rnacells)
CCLE_rna <- summarizeMolecularProfiles(CCLE, mDataType = "Kallisto_0.46.1.isoforms.counts", cell.lines = intersected_rnacells)
GDSC_rna <- summarizeMolecularProfiles(GDSC, mDataType = "Kallisto_0.46.1.isoforms.counts", cell.lines = intersected_rnacells)

#load("../results/Kallisto_isoforms_counts.RData")

##keep protein-coding transcripts only to ensure consistency (N=84,107). Protein transcripts obtained 
transcript_gene <- GDSC_rna@elementMetadata[,c("gene_id", "transcript_id")] #save corresponding gene_id for gene expression quantification 
protein_transcripts <-gCSI_rna@elementMetadata$transcript_id[which(gCSI_rna@elementMetadata$transcript_type == "protein_coding")]

expr_gcsi_i <- gCSI_rna@assays@data$expr[protein_transcripts,]
expr_ccle_i <- CCLE_rna@assays@data$expr[protein_transcripts,]
expr_gdsc_i <- GDSC_rna@assays@data$expr[protein_transcripts,]

##keep only transcripts that do not have 0 expression across all cell lines (N=80,786)
expr_gcsi_i <- expr_gcsi_i[apply(expr_gcsi_i[,-1], 1, function(x) !all(x==0)),]
expr_ccle_i <- expr_ccle_i[apply(expr_ccle_i[,-1], 1, function(x) !all(x==0)),]
expr_gdsc_i <- expr_gdsc_i[apply(expr_gdsc_i[,-1], 1, function(x) !all(x==0)),]

##only keep transcripts found in all datasets now after filitering for 0 expression above
intersected_t <- Reduce(intersect, list(rownames(expr_gcsi_i), rownames(expr_ccle_i), rownames(expr_gdsc_i)))
expr_gcsi_i <- expr_gcsi_i[intersected_t,]
expr_ccle_i <- expr_ccle_i[intersected_t,]
expr_gdsc_i <- expr_gdsc_i[intersected_t,]

#remove housekeeping transcripts (defined by HRT Atlas)
transcripts_noid <- gsub("\\..*","",intersected_t)
housekeeping <- read.csv("/data/transcript_stability/Housekeeping_TranscriptsHuman.csv", sep = ";")

expr_gcsi_i <- expr_gcsi_i[which(!transcripts_noid %in%  housekeeping$Ensembl),]
expr_ccle_i <- expr_ccle_i[which(!transcripts_noid %in%  housekeeping$Ensembl),]
expr_gdsc_i <- expr_gdsc_i[which(!transcripts_noid %in%  housekeeping$Ensembl),]

save(expr_gcsi_i, expr_ccle_i, expr_gdsc_i, file = "../results/isoform_expression.RData")


######################################
##### Get gene transcript counts #####
######################################

#summarize gene expression across replicates for each pSet
gCSI_rna <- summarizeMolecularProfiles(gCSI, mDataType = "Kallisto_0.46.1.rnaseq.counts", cell.lines = intersected_rnacells)
CCLE_rna <- summarizeMolecularProfiles(CCLE, mDataType = "Kallisto_0.46.1.rnaseq.counts", cell.lines = intersected_rnacells)
GDSC_rna <- summarizeMolecularProfiles(GDSC, mDataType = "Kallisto_0.46.1.rnaseq.counts", cell.lines = intersected_rnacells)

#load("../results/Kallisto_rnaseq_counts.RData")

##keep protein-coding genes that correspond to transcripts from transcript_stability.R only to ensure consistency (N=19,957). Protein genes obtained 
protein_genes <-GDSC_rna@elementMetadata$gene_id[which(GDSC_rna@elementMetadata$gene_type == "protein_coding" & GDSC_rna@elementMetadata$gene_id %in% transcript_gene$gene_id)]

expr_gcsi_p <- gCSI_rna@assays@data$expr[protein_genes,]
expr_ccle_p <- CCLE_rna@assays@data$expr[protein_genes,]
expr_gdsc_p <- GDSC_rna@assays@data$expr[protein_genes,]


##keep only genes that do not have 0 expression across all cell lines (N=80,786)
expr_gcsi_p <- expr_gcsi_p[apply(expr_gcsi_p[,-1], 1, function(x) !all(x==0)),]
expr_ccle_p <- expr_ccle_p[apply(expr_ccle_p[,-1], 1, function(x) !all(x==0)),]
expr_gdsc_p <- expr_gdsc_p[apply(expr_gdsc_p[,-1], 1, function(x) !all(x==0)),]

##only keep genes found in all datasets now after filitering for 0 expression above (N = 19957)
intersected_g <- Reduce(intersect, list(rownames(expr_gcsi_p), rownames(expr_ccle_p), rownames(expr_gdsc_p)))
expr_gcsi_p <- expr_gcsi_p[intersected_g,]
expr_ccle_p <- expr_ccle_p[intersected_g,]
expr_gdsc_p <- expr_gdsc_p[intersected_g,]

save(expr_gcsi_p, expr_ccle_p, expr_gdsc_p, file = "../results/gene_expression.RData")


