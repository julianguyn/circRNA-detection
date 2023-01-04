##compares predictiveness of isoforms for a known biomarker (gene) among ranked stability indicies (are the most "stable" isoforms, the most predictive?)

suppressMessages(library(PharmacoGx))
suppressMessages(library(survcomp))
suppressMessages(library(gimme))
suppressMessages(library(dplyr))
suppressMessages(library(forestplot))
suppressMessages(library(reshape2))

#load biological replicates
load("../results/inter.RData")

#read in PSets
gCSI <- readRDS("/data/PSets/gCSI.rds")
CCLE <- readRDS("/data/PSets/CCLE.rds")
GDSC <- readRDS("/data/PSets/GDSC2-8.2.rds")

#update PSets for latest PharmacoGx update
gCSI <- updateObject(gCSI)
CCLE <- updateObject(CCLE)
GDSC <- updateObject(GDSC)

#canonical transcripts
canonical <- read.csv("/data/biomarker/canonical_transcripts", sep="\t")
canonical$t_name <- gCSI@molecularProfiles$Kallisto_0.46.1.isoforms@elementMetadata$transcript_name[match(canonical$transcript, gCSI@molecularProfiles$Kallisto_0.46.1.isoforms@elementMetadata$transcript_id)]

#modify annotations so that the molecular profiles can be read in the permutation function
S4Vectors::metadata(gCSI@molecularProfiles[["Kallisto_0.46.1.isoforms.counts"]])$annotation <- "isoform"
S4Vectors::metadata(CCLE@molecularProfiles[["Kallisto_0.46.1.isoforms.counts"]])$annotation <- "isoform"
S4Vectors::metadata(GDSC@molecularProfiles[["Kallisto_0.46.1.isoforms.counts"]])$annotation <- "isoform"

#remove tissue id from cell to prevent Simpsons paradox from affecting pearson computation
sampleInfo(gCSI)$tissueid <- NA
sampleInfo(CCLE)$tissueid <- NA
sampleInfo(GDSC)$tissueid <- NA

#function that creates data.frame of every single combination between two vectors (e.g. drug and gene combinations)
expand.grid.unique <- function(x, y, include.equals=FALSE)
{
  x <- unique(x)
  y <- unique(y)
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}

#read in known biomarkers (curated by BHKLAB), and filter for expression biomarkers only with drugs that exist across gCSI,CCLE,GDSC
#biomarkers selected are: 1) ERBB2 + Lapatinib; 2) ALK + Crizotinib; 3) PHB + Paclitaxel; 4) ESR2 + Erlotinib; 5) EGFR & Erlotinib

intersected_drugs <- Reduce(intersect, list(colnames(gCSI@treatmentResponse$n), colnames(CCLE@treatmentResponse$n), colnames(GDSC@treatmentResponse$n)))
genes <- c("ERBB2","ALK","PHB","ESR2","EGFR")
drugs <- c("Lapatinib", "Crizotinib", "Paclitaxel", "Erlotinib", "Erlotinib")
biomarkers <- data.frame(matrix(data = NA, nrow = length(genes), ncol = 2))
colnames(biomarkers) <- c("compound","gene")
biomarkers$compound <- drugs
biomarkers$gene <- genes
gene_ids <- gCSI@molecularProfiles$Kallisto_0.46.1.rnaseq@elementMetadata$gene_id[match(biomarkers$gene, gCSI@molecularProfiles$Kallisto_0.46.1.rnaseq@elementMetadata$gene_name)]
biomarkers$gene_ids <- gene_ids

####################
####COMPUTE CI #####
####################


computeCI <- function(PSet, mData, features, cells, drugs){
  
  rnaseq <- summarizeMolecularProfiles(PSet, mDataType = mData)
  rnaseq <- as.data.frame(rnaseq@assays@data$exprs)
  
  sensitivity_data <- as.data.frame(summarizeSensitivityProfiles(PSet, sensitivity.measure = "aac_recomputed",  fill.missing = F))
  intersected_rnacells <- cells
  commonSamples <- intersected_rnacells
  commonSamples <- intersected_rnacells[which(intersected_rnacells %in% colnames(sensitivity_data))]
  combinations <- as.data.frame(matrix(data = NA, nrow = length(drugs), ncol = 3))
  if (!mData == "Kallisto_0.46.1.isoforms.counts"){
    gene_ids <- gCSI@molecularProfiles$Kallisto_0.46.1.rnaseq@elementMetadata$gene_id[match(features, gCSI@molecularProfiles$Kallisto_0.46.1.rnaseq@elementMetadata$gene_name)]
    colnames(combinations) <- c("gene","drug", "gene_id")
    combinations$gene <- features
    combinations$drug <- drugs
    combinations$gene_id <- gene_ids
    combinations$ci <- NA
    combinations$pvalue <- NA
    combinations$se <- NA
    combinations$upper <- NA
    combinations$lower <- NA
  }else{
    transcript_name <- gCSI@molecularProfiles$Kallisto_0.46.1.isoforms.counts@elementMetadata$transcript_name[match(features, gCSI@molecularProfiles$Kallisto_0.46.1.isoforms.counts@elementMetadata$transcript_id)]
    colnames(combinations) <- c("transcript","drug", "transcript_id")
    combinations$transcript <- transcript_name
    combinations$drug <- drugs
    combinations$transcript_id <- features
    combinations$ci <- NA
    combinations$pvalue <- NA
    combinations$se <- NA
    combinations$upper <- NA
    combinations$lower <- NA
  }
  
  for (i in 1:nrow(combinations)){
    print(paste0(i, " out of ", nrow(combinations), " complete"))
    tt <- sensitivity_data[combinations[,2][i], commonSamples]
    ci <- survcomp::concordance.index(as.numeric(tt), surv.time = as.numeric(unlist(-rnaseq[combinations[,3][i], commonSamples])), 
                                      surv.event = rep(1,length(sensitivity_data[commonSamples])), 
                                      outx = TRUE, method="noether", na.rm = TRUE)
    
    combinations$pvalue[i] <- ci$p.value
    combinations$ci[i] <- ci$c.index
    combinations$se[i] <- ci$se
    combinations$upper[i] <- ci$upper
    combinations$lower[i] <- ci$lower
  }
  
  return(combinations)
  
}

gCSI_CI_gene <- computeCI(PSet= gCSI, 
                          mData = "Kallisto_0.46.1.rnaseq", 
                          cells = intersected_rnacells,
                          features = biomarkers$gene,
                          drugs = biomarkers$compound)

CCLE_CI_gene <- computeCI(PSet= CCLE, 
                          mData = "Kallisto_0.46.1.rnaseq", 
                          cells = intersected_rnacells,
                          features = biomarkers$gene,
                          drugs = biomarkers$compound)

GDSC_CI_gene <- computeCI(PSet= GDSC, 
                          mData = "Kallisto_0.46.1.rnaseq", 
                          cells = intersected_rnacells,
                          features = biomarkers$gene,
                          drugs = biomarkers$compound)


####################################
#####COMPUTE CI FOR TRANSCRIPTS#####
####################################

##import transcript stability data generated by transcript_stability.R
transcript_stability <- read.csv("../results/transcript_stability.csv")
biomarkers_transcript <- list()
for (i in biomarkers$gene){
  
  t_id <- gCSI@molecularProfiles$Kallisto_0.46.1.isoforms@elementMetadata$transcript_id[which(gCSI@molecularProfiles$Kallisto_0.46.1.isoforms@elementMetadata$gene_name == i)]
  t_names_matched <- transcript_stability$transcript_id[which(transcript_stability$transcript_id %in% t_id)]
  d <- biomarkers$compound[which(biomarkers$gene == i)]
  df <- as.data.frame(expand.grid.unique(t_names_matched, d,include.equals = TRUE))
  biomarkers_transcript[[i]] <- df
  
}

biomarkers_transcript <- melt(biomarkers_transcript)
colnames(biomarkers_transcript) <- c("transcript_id", "compound", "gene")


gCSI_CI_t <- computeCI(PSet= gCSI, 
                       mData = "Kallisto_0.46.1.isoforms.counts", 
                       cells = intersected_rnacells,
                       features = biomarkers_transcript$transcript_id,
                       drugs = biomarkers_transcript$compound)

CCLE_CI_t <- computeCI(PSet= CCLE, 
                       mData = "Kallisto_0.46.1.isoforms.counts", 
                       cells = intersected_rnacells,
                       features = biomarkers_transcript$transcript_id,
                       drugs = biomarkers_transcript$compound)

GDSC_CI_t <- computeCI(PSet= GDSC, 
                       mData = "Kallisto_0.46.1.isoforms.counts", 
                       cells = intersected_rnacells,
                       features = biomarkers_transcript$transcript_id,
                       drugs = biomarkers_transcript$compound)


save(gCSI_CI_gene, CCLE_CI_gene, GDSC_CI_gene, file="../results/gene_ci.RData")
save(gCSI_CI_t, CCLE_CI_t, GDSC_CI_t, file="../results/transcript_ci.RData")

#########################################
####COMPUTE PEARSON (w/Permutations)#####
#########################################

computePearson <- function(PSet, mData, features, cells, drugs){
  
  combinations <- as.data.frame(matrix(data = NA, nrow = length(drugs), ncol = 3))
  if (!mData == "Kallisto_0.46.1.isoforms.counts"){
    gene_ids <- gCSI@molecularProfiles$Kallisto_0.46.1.rnaseq@elementMetadata$gene_id[match(features, gCSI@molecularProfiles$Kallisto_0.46.1.rnaseq@elementMetadata$gene_name)]
    colnames(combinations) <- c("gene","drug", "gene_id")
    combinations$gene <- features
    combinations$drug <- drugs
    combinations$gene_id <- gene_ids
    combinations$pearson <- NA
    combinations$pvalue <- NA
    combinations$upper <- NA
    combinations$lower <- NA
    combinations$se <- NA
  }else{
    transcript_name <- gCSI@molecularProfiles$Kallisto_0.46.1.isoforms.counts@elementMetadata$transcript_name[match(features, gCSI@molecularProfiles$Kallisto_0.46.1.isoforms.counts@elementMetadata$transcript_id)]
    colnames(combinations) <- c("transcript","drug", "transcript_id")
    combinations$transcript <- transcript_name
    combinations$drug <- drugs
    combinations$transcript_id <- features
    combinations$pearson <- NA
    combinations$pvalue <- NA
    combinations$upper <- NA
    combinations$lower <- NA
    combinations$se <- NA
  }
  
  for (i in 1:nrow(combinations)){
    print(paste0(i, " out of ", nrow(combinations), " complete"))
    
    sig <- drugSensitivitySig(PSet, mData, drugs=combinations[i,2], 
                              features=combinations[i,3], sensitivity.measure = "aac_recomputed", modeling.method="pearson", 
                              inference.method="resampling", cells=cells, nthread=4, parallel.on = "drug")
    
    combinations$pvalue[i] <- sig[,,"pvalue"]
    combinations$pearson[i] <- sig[,,"estimate"]
    combinations$upper[i] <- sig[,,"upper"]
    combinations$lower[i] <- sig[,,"lower"]
    combinations$se[i] <- (combinations$upper[i] - combinations$lower[i])/3.92
  }
  
  return(combinations)
  
}

#permutation function not yet publically availble in PharmacoGx.

#gCSI_PEAR_gene <- computePearson(PSet= gCSI, 
#mData = "Kallisto_0.46.1.rnaseq", 
#cells = intersected_rnacells,
#features = biomarkers$gene,
#drugs = biomarkers$compound)


#CCLE_PEAR_gene <- computePearson(PSet= CCLE, 
#mData = "Kallisto_0.46.1.rnaseq", 
#cells = intersected_rnacells,
#features = biomarkers$gene,
#drugs = biomarkers$compound)


#GDSC_PEAR_gene <- computePearson(PSet= GDSC, 
#mData = "Kallisto_0.46.1.rnaseq", 
#cells = intersected_rnacells,
#features = biomarkers$gene,
#drugs = biomarkers$compound)



#gCSI_PEAR_t <- computePearson(PSet= gCSI, 
#mData = "Kallisto_0.46.1.isoforms.counts", 
#cells = intersected_rnacells,
#features = biomarkers_transcript$transcript_id,
#drugs = biomarkers_transcript$compound)

#CCLE_PEAR_t <- computePearson(PSet= CCLE, 
#mData = "Kallisto_0.46.1.isoforms.counts", 
#cells = intersected_rnacells,
#features = biomarkers_transcript$transcript_id,
#drugs = biomarkers_transcript$compound)

#GDSC_PEAR_t <- computePearson(PSet= GDSC, 
#mData = "Kallisto_0.46.1.isoforms.counts", 
#cells = intersected_rnacells,
#features = biomarkers_transcript$transcript_id,
# drugs = biomarkers_transcript$compound)


load("/data/biomarker/transcript_pearson.RData")
load("/data/biomarker/gene_pearson.RData")

###############################
##Compute mean gene stability##
###############################

CCLE_rna <- summarizeMolecularProfiles(CCLE, mData="Kallisto_0.46.1.rnaseq", intersected_rnacells)
gCSI_rna <- summarizeMolecularProfiles(gCSI, mData="Kallisto_0.46.1.rnaseq", intersected_rnacells)
GDSC_rna <- summarizeMolecularProfiles(GDSC, mData="Kallisto_0.46.1.rnaseq", intersected_rnacells)

erbb2 <- "ENSG00000141736.13"
alk <- "ENSG00000171094.18"
egfr <- "ENSG00000146648.18"
phb <- "ENSG00000167085.11"
esr2 <- "ENSG00000140009.18"

gene_ids <- c(erbb2, alk, egfr, phb, esr2)

gene_stability <- data.frame(matrix(ncol = 1, nrow=5))
colnames(gene_stability) <- c("mean_stability")
rownames(gene_stability) <- gene_ids

for (i in gene_ids){
  
  e1 <- CCLE_rna@assays@data$exprs[i, ]
  e2 <- gCSI_rna@assays@data$exprs[i, ]
  e3 <- GDSC_rna@assays@data$exprs[i, ]
  
  gcsi_ccle <- cor(e2,e1, method="spearman")
  gcsi_gdsc <- cor(e2,e3, method="spearman")
  gdsc_ccle <- cor(e3,e1, method="spearman")
  
  gene_stability[i, "mean_stability"] <- mean(gcsi_ccle, gcsi_gdsc, gdsc_ccle)
  
}

rownames(gene_stability) <- c("ERBB2", "ALK","EGFR", "PHB","ESR2")

#####################################
##Compute mean transcript stability##
#####################################

gCSI_CI_t$mean_stability <- as.numeric(rowMeans(transcript_stability[match(gCSI_CI_t$transcript_id, transcript_stability$transcript_id),
                                                                     c("gcsi_ccle_spearman", "gcsi_gdsc_spearman","gdsc_ccle_spearman")]))

CCLE_CI_t$mean_stability <- as.numeric(rowMeans(transcript_stability[match(CCLE_CI_t$transcript_id, transcript_stability$transcript_id),
                                                                     c("gcsi_ccle_spearman", "gcsi_gdsc_spearman","gdsc_ccle_spearman")]))

GDSC_CI_t$mean_stability <- as.numeric(rowMeans(transcript_stability[match(GDSC_CI_t$transcript_id, transcript_stability$transcript_id),
                                                                     c("gcsi_ccle_spearman", "gcsi_gdsc_spearman","gdsc_ccle_spearman")]))



gCSI_PEAR_t$mean_stability <- as.numeric(rowMeans(transcript_stability[match(gCSI_PEAR_t$transcript_id, transcript_stability$transcript_id),
                                                                       c("gcsi_ccle_spearman", "gcsi_gdsc_spearman","gdsc_ccle_spearman")]))

CCLE_PEAR_t$mean_stability <- as.numeric(rowMeans(transcript_stability[match(CCLE_PEAR_t$transcript_id, transcript_stability$transcript_id),
                                                                       c("gcsi_ccle_spearman", "gcsi_gdsc_spearman","gdsc_ccle_spearman")]))

GDSC_PEAR_t$mean_stability <- as.numeric(rowMeans(transcript_stability[match(GDSC_PEAR_t$transcript_id, transcript_stability$transcript_id),
                                                                       c("gcsi_ccle_spearman", "gcsi_gdsc_spearman","gdsc_ccle_spearman")]))


dir.create("../results/CI_forest")
dir.create("../results/Pearson_forest")
dir.create("../results/Meta_forest")

##################################################
##PLOTS 1: GENE(ISOFORM)/DRUG PER DATASET -> CI ##
##################################################

dataset = c("gCSI","CCLE","GDSC")

for (d in dataset){
  
  if (d == "gCSI"){
    d_gene = gCSI_CI_gene 
    d_transcript = gCSI_CI_t 
  } else if (d == "CCLE") {
    d_gene = CCLE_CI_gene 
    d_transcript = CCLE_CI_t 
  } else {
    d_gene = GDSC_CI_gene 
    d_transcript = GDSC_CI_t 
    
  }
  
  for (i in 1:nrow(d_gene)){
    gene <- d_gene$gene[i]
    drug <- d_gene$drug[i]
    
    transcripts <- d_transcript[grep(gene, d_transcript$transcript),]
    transcripts <- transcripts[order(transcripts$mean_stability, decreasing = TRUE),]
    
    gene_ci_value <- d_gene[which(d_gene$gene == gene),"ci"]
    gene_ci_pvalue <- d_gene[which(d_gene$gene == gene),"pvalue"]
    gene_ci_lower <-   d_gene[which(d_gene$gene == gene),"lower"]
    gene_ci_upper <-   d_gene[which(d_gene$gene == gene),"upper"]
    
    ll1 <- nrow(transcripts) + 2
    c_indices <- structure(
      list(
        mean  = c(NA, gene_ci_value, transcripts$ci),
        lower = c(NA, gene_ci_lower, transcripts$lower),
        upper = c(NA, gene_ci_upper, transcripts$upper)
      ),
      .Names = c("C-index    ", "lower", "upper"),
      row.names = c(NA, -(ll1)), 
      class = "data.frame"
    )
    
    
    c_tabletext <- cbind(
      c("Gene/Transcript", paste0(gene, " (Gene)"), c(transcripts$transcript)),
      c("N",rep(144,length(transcripts$transcript)+1)),
      
      c("C-index", formatC(gene_ci_value, digits = 2),
        formatC(transcripts$ci, digits = 2)),
      
      c("P-value", formatC(gene_ci_pvalue, format = "e",digits = 2),
        formatC(transcripts$pvalue, format = "e", digits = 2)),
      
      
      c("Mean stability", formatC(gene_stability[gene,"mean_stability"], digits = 2), formatC(transcripts$mean_stability, digits = 2))
    )
    
    
    fn <- local({
      i = 0
      
      b_clrs =  c(c("darkgreen"),c(rep("red", nrow(transcripts))))
      l_clrs =  c(c("darkgreen"),c(rep("red", nrow(transcripts))))
      function(..., clr.line, clr.marker){
        i <<- i + 1
        fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
        #fpDrawSummaryCI(...,col=s_clrs[i])
      }
    })
    
    
    
    t <- c_tabletext[3:nrow(c_tabletext),1]
    t_matched <- which(t %in% canonical$t_name)
    if (length(t_matched) == 0){
      c_t <- fpColors(text="black")
    } else {
      ll <- nrow(c_tabletext) - 2
      cc <- rep("black", ll)
      cc[t_matched] <- "blue"
      c_t <- fpColors(text=c("black","black",cc))
    }
    
    c_index_matrix <- list()
    
    c_index_matrix <- c(list(c_tabletext), c_index_matrix)
    fileName = paste0("../results/CI_forest/",d,"_CI_",gene,"_",drug,".pdf")
    pdf(fileName, width=10 , height=8, onefile=FALSE)
    
    ll <- nrow(c_tabletext) - 1
    forestplot(c_tabletext, c_indices, new_page = TRUE, boxsize = 0.3, is.summary=c(T,rep(F,ll)), xlab="Concordance Index", 
               title=paste0(gene,"_",drug," (",d,")"), zero=c(.49, .51),hrzl_lines=list("2"=gpar(lty=2, columns=1:5, col = "#000044")),
               txt_gp=fpTxtGp(label=gpar(fontfamily = "", cex = 1.1, fontface=2),
                              ticks=gpar(fontfamily = "", cex=.5, fontface=1),
                              xlab=gpar(fontfamily = "", cex=0.9, fontface=2), 
                              legend=gpar(fontfamily = "", cex = 0.5, fontface=1)),
               fn.ci_norm = fn,
               col = c_t,
               xticks= c(.3, .4, .5, .6, .7, .8)
    )
    
    dev.off()
    
  }
  
}


#######################################################
##PLOTS 2: GENE(ISOFORM)/DRUG PER DATASET -> PEARSON ##
#######################################################


dataset = c("gCSI","CCLE","GDSC")

for (d in dataset){
  
  if (d == "gCSI"){
    d_gene = gCSI_PEAR_gene 
    d_transcript = gCSI_PEAR_t 
  } else if (d == "CCLE") {
    d_gene = CCLE_PEAR_gene
    d_transcript = CCLE_PEAR_t 
  } else {
    d_gene = GDSC_PEAR_gene
    d_transcript = GDSC_PEAR_t
    
  }
  
  for (i in 1:nrow(d_gene)){
    gene <- d_gene$gene[i]
    drug <- d_gene$drug[i]
    
    transcripts <- d_transcript[grep(gene, d_transcript$transcript),]
    transcripts <- transcripts[order(transcripts$mean_stability, decreasing = TRUE),]
    
    gene_pear_value <- d_gene[which(d_gene$gene == gene),"pearson"]
    gene_pear_pvalue <- d_gene[which(d_gene$gene == gene),"pvalue"]
    gene_pear_lower <-   d_gene[which(d_gene$gene == gene),"lower"]
    gene_pear_upper <-   d_gene[which(d_gene$gene == gene),"upper"]
    
    ll1 <- nrow(transcripts) + 2
    c_indices <- structure(
      list(
        mean  = c(NA, gene_pear_value, transcripts$pearson),
        lower = c(NA, gene_pear_lower, transcripts$lower),
        upper = c(NA, gene_pear_upper, transcripts$upper)
      ),
      .Names = c("C-index    ", "lower", "upper"),
      row.names = c(NA, -(ll1)), 
      class = "data.frame"
    )
    
    
    c_tabletext <- cbind(
      c("Gene/Transcript", paste0(gene, " (Gene)"), c(transcripts$transcript)),
      c("N",rep(144,length(transcripts$transcript)+1)),
      
      c("Pearson", formatC(gene_pear_value, digits = 2),
        formatC(transcripts$pearson, digits = 2)),
      
      c("P-value", formatC(gene_pear_pvalue, format = "e", digits = 2),
        formatC(transcripts$pvalue, format = "e", digits = 2)),
      
      
      c("Mean stability", formatC(gene_stability[gene,"mean_stability"], digits = 2), formatC(transcripts$mean_stability, digits = 2))
    )
    
    
    fn <- local({
      i = 0
      
      b_clrs =  c(c("darkgreen"),c(rep("red", nrow(transcripts))))
      l_clrs =  c(c("darkgreen"),c(rep("red", nrow(transcripts))))
      function(..., clr.line, clr.marker){
        i <<- i + 1
        fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
        #fpDrawSummaryCI(...,col=s_clrs[i])
      }
    })
    
    
    
    t <- c_tabletext[3:nrow(c_tabletext),1]
    t_matched <- which(t %in% canonical$t_name)
    if (length(t_matched) == 0){
      c_t <- fpColors(text="black")
    } else {
      ll <- nrow(c_tabletext) - 2
      cc <- rep("black", ll)
      cc[t_matched] <- "blue"
      c_t <- fpColors(text=c("black","black",cc))
    }
    
    c_index_matrix <- list()
    
    c_index_matrix <- c(list(c_tabletext), c_index_matrix)
    fileName = paste0("../results/Pearson_forest/",d,"_Pearson_",gene,"_",drug,".pdf")
    pdf(fileName, width=10 , height=8, onefile=FALSE)
    
    ll <- nrow(c_tabletext) - 1
    forestplot(c_tabletext, c_indices, new_page = TRUE, boxsize = 0.3, is.summary=c(T,rep(F,ll)), xlab="Pearson Coefficient", 
               title=paste0(gene,"_",drug," (",d,")"), zero=c(-.04, .04),hrzl_lines=list("2"=gpar(lty=2, columns=1:5, col = "#000044")),
               txt_gp=fpTxtGp(label=gpar(fontfamily = "", cex = 1.1, fontface=2),
                              ticks=gpar(fontfamily = "", cex=.5, fontface=1),
                              xlab=gpar(fontfamily = "", cex=0.9, fontface=2), 
                              legend=gpar(fontfamily = "", cex = 0.5, fontface=1)),
               fn.ci_norm = fn,
               col = c_t,
               xticks= c(-1, -.5, 0, .5, 1)
    )
    
    dev.off()
    
  }
  
}


###################################################
##PLOT 3: GENE(ISOFORM)/DRUG PER DATASET -> META ##
###################################################

for (i in 1:length(genes)){
  gene <- genes[i]
  drug <- drugs[i]
  
  gene_CI <- cbind(gcsi=gCSI_CI_gene, ccle=CCLE_CI_gene, gdsc=GDSC_CI_gene)
  colnames(gene_CI)[1:2] <- c("gene","drug")
  
  gene_permut <- cbind(gcsi=gCSI_PEAR_gene, ccle=CCLE_PEAR_gene, gdsc=GDSC_PEAR_gene)
  colnames(gene_permut)[1:2] <- c("gene","drug")
  
  gene.CI <- gene_CI[which(gene_CI$gene == gene & gene_CI$drug == drug),]
  gene_perm <- gene_permut[which(gene_permut$gene == gene & gene_permut$drug == drug),]
  
  isoform_CI <- cbind(gcsi=gCSI_CI_t, ccle=CCLE_CI_t, gdsc=GDSC_CI_t)
  isoform_CI$gene <- gsub("\\-.*","",isoform_CI$gcsi.transcript)
  colnames(isoform_CI)[2] <- "drug"
  
  isoform_perm <- cbind(gcsi=gCSI_PEAR_t, ccle=CCLE_PEAR_t, gdsc=GDSC_PEAR_t)
  isoform_perm$gene <- gsub("\\-.*","",isoform_perm$gcsi.transcript)
  colnames(isoform_perm)[2] <- "drug"
  
  
  transcript_CI <- isoform_CI[which(isoform_CI$gene == gene & isoform_CI$drug == drug),]  
  transcript_perm <- isoform_perm[which(isoform_perm$gene == gene & isoform_perm$drug == drug),]
  
  
  #order transcripts based on stability (mean stability)
  transcript_CI[,c("gdsc.mean_stability","ccle.mean_stability")] <- NULL
  transcript_perm[,c("gdsc.mean_stability","ccle.mean_stability")] <- NULL
  
  transcript_CI <- transcript_CI[order(transcript_CI$gcsi.mean_stability, decreasing = TRUE),] 
  transcript_perm <- transcript_perm[order(transcript_perm$gcsi.mean_stability, decreasing = TRUE),] 
  
  
  plot_combined_meta <- data.frame(matrix(nrow = nrow(transcript_perm), ncol=9))
  colnames(plot_combined_meta) <- c("gene","ci_estimate", "ci_se", "ci_upper","ci_lower", "perm_estimate","perm_se","perm_upper","perm_lower")
  plot_combined_meta$gene <- transcript_perm$gene
  
  #combined gene CI
  combined_ci <- combine.est(
    c(
      gene.CI$gcsi.ci, gene.CI$ccle.ci, gene.CI$gdsc.ci
    ),
    c(
      gene.CI$gcsi.se, gene.CI$ccle.se, gene.CI$gdsc.se
    ),na.rm = TRUE,hetero = TRUE)
  
  #combined gene pearson correlation from permutations
  combined_pear <- combine.est(
    c(
      gene_perm$gcsi.pearson, gene_perm$ccle.pearson, gene_perm$gdsc.pearson
    ),
    c(
      gene_perm$gcsi.se, gene_perm$ccle.se, gene_perm$gdsc.se
    ),na.rm = TRUE,hetero = TRUE)
  
  #loop for combined isoform pearson correlation from permutations
  for (i in 1:nrow(transcript_perm)){
    plot_combined_meta$ci_estimate[i] <- combine.est(
      c(
        transcript_CI$gcsi.ci[i], transcript_CI$ccle.ci[i], transcript_CI$gdsc.ci[i]
      ),
      c(
        transcript_CI$gcsi.se[i], transcript_CI$ccle.se[i], transcript_CI$gdsc.se[i]
      ),na.rm = TRUE,hetero = TRUE)$estimate
    
    plot_combined_meta$ci_se[i] <- combine.est(
      c(
        transcript_CI$gcsi.ci[i], transcript_CI$ccle.ci[i], transcript_CI$gdsc.ci[i]
      ),
      c(
        transcript_CI$gcsi.se[i], transcript_CI$ccle.se[i], transcript_CI$gdsc.se[i]
      ),na.rm = TRUE,hetero = TRUE)$se
    
    
    plot_combined_meta$perm_estimate[i] <- combine.est(
      c(
        transcript_perm$gcsi.pearson[i], transcript_perm$ccle.pearson[i], transcript_perm$gdsc.pearson[i]
      ),
      c(
        transcript_perm$gcsi.se[i], transcript_perm$ccle.se[i], transcript_perm$gdsc.se[i]
      ),na.rm = TRUE,hetero = TRUE)$estimate
    
    
    plot_combined_meta$perm_se[i] <- combine.est(
      c(
        transcript_perm$gcsi.pearson[i], transcript_perm$ccle.pearson[i], transcript_perm$gdsc.pearson[i]
      ),
      c(
        transcript_perm$gcsi.se[i], transcript_perm$ccle.se[i], transcript_perm$gdsc.se[i]
      ),na.rm = TRUE,hetero = TRUE)$se
    
    plot_combined_meta$ci_lower[i] <- plot_combined_meta$ci_estimate[i] + qnorm(0.025, lower.tail=TRUE) *  plot_combined_meta$ci_se[i]
    plot_combined_meta$ci_upper[i] <- plot_combined_meta$ci_estimate[i] + qnorm(0.025, lower.tail=FALSE) *  plot_combined_meta$ci_se[i]
    
    plot_combined_meta$perm_lower[i] <- plot_combined_meta$perm_estimate[i] + qnorm(0.025, lower.tail=TRUE) *  plot_combined_meta$perm_se[i]
    plot_combined_meta$perm_upper[i] <- plot_combined_meta$perm_estimate[i] + qnorm(0.025, lower.tail=FALSE) *  plot_combined_meta$perm_se[i]
    
  }
  
  
  combined_ci_lower <- combined_ci$estimate + qnorm(0.025, lower.tail=TRUE) *  combined_ci$se
  combined_ci_upper <- combined_ci$estimate + qnorm(0.025, lower.tail=FALSE) *  combined_ci$se
  
  combined_pear_lower <- combined_pear$estimate + qnorm(0.025, lower.tail=TRUE) *  combined_pear$se
  combined_pear_upper <- combined_pear$estimate + qnorm(0.025, lower.tail=FALSE) *  combined_pear$se
  
  #combined_ci_p <- pnorm((combined_ci$estimate - 0.5)/combined_ci$se, lower.tail = combined_ci$estimate < 0.5) * 2
  
  rows_l <- length(c(NA, combined_ci$estimate,
                     plot_combined_meta$ci_estimate))
  
  c_indices <- structure(
    list(
      mean  = c(NA, combined_ci$estimate,
                plot_combined_meta$ci_estimate),
      
      lower = c(NA, combined_ci_lower,
                plot_combined_meta$ci_lower),
      
      upper = c(NA, combined_ci_upper,
                plot_combined_meta$ci_upper)
    ),
    .Names = c("C-index    ", "lower", "upper"),
    row.names = c(NA, -(rows_l)), 
    class = "data.frame"
  )
  
  
  c_tabletext <- cbind(
    c("Gene/Transcript", paste0(gene, " (Gene)"), c(transcript_CI$gcsi.transcript)),
    c("N",rep(144,rows_l-1)),
    
    c("C-index", formatC(combined_ci$estimate, digits = 2, format = "f"),
      formatC(plot_combined_meta$ci_estimate, digits = 2, format = "f")),
    
    c("Pearson", formatC(combined_pear$estimate, digits = 2, format = "f"),
      formatC(plot_combined_meta$perm_estimate, digits = 2, format = "f")),
    
    
    c("Mean stability", formatC(gene_stability[gene,"mean_stability"], digits = 2, format = "f"), formatC(transcript_perm$gcsi.mean_stability, digits = 2, format = "f"))
  )
  
  
  fn <- local({
    i = 0
    
    b_clrs =  c(c("darkgreen",rep("red", nrow(transcript_CI))))
    l_clrs =  c(c("darkgreen",rep("red", nrow(transcript_CI))))
    function(..., clr.line, clr.marker){
      i <<- i + 1
      fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
      #fpDrawSummaryCI(...,col=s_clrs[i])
    }
  })
  
  fn1 <- local({
    i = 0
    
    s_clrs = c(c("darkgreen",rep("red", nrow(transcript_CI))))
    function(..., col){
      i <<- i + 1
      fpDrawSummaryCI(...,col=s_clrs[i])
    }
  })
  
  
  t <- c_tabletext[3:nrow(c_tabletext),1]
  t_matched <- which(t %in% canonical$t_name)
  if (length(t_matched) == 0){
    c_t <- fpColors(text="black")
  } else {
    ll <- nrow(c_tabletext) - 2
    cc <- rep("black", ll)
    cc[t_matched] <- "blue"
    c_t <- fpColors(text=c("black","black",cc))
  }
  
  c_index_matrix <- list()
  
  c_index_matrix <- c(list(c_tabletext), c_index_matrix)
  
  fileName = paste0("../results/Meta_forest/",gene,"_",drug,".pdf")
  pdf(fileName, width=10 , height=8, onefile=FALSE)
  
  forestplot(c_tabletext, c_indices, new_page = TRUE, boxsize = 0.3, is.summary=c(T), xlab="Concordance Index", 
             title=paste0(gene,"_",drug), zero=c(.49, .51),hrzl_lines=list("2"=gpar(lty=2, columns=1:5, col = "#000044")),
             txt_gp=fpTxtGp(label=gpar(fontfamily = "", cex = 1.1, fontface=2),
                            ticks=gpar(fontfamily = "", cex=1.1, fontface=1),
                            xlab=gpar(fontfamily = "", cex=0.9, fontface=2), 
                            legend=gpar(fontfamily = "", cex = 0.5, fontface=1)),
             fn.ci_sum = fn1,
             fn.ci_norm = fn,
             col = c_t,
             xticks= c(.3, .4, .5, .6, .7, .8)
  )
  
  dev.off()
  
}
