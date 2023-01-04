suppressMessages(library(genefu))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))



transcript_stability <- read.csv("../results/transcript_stability.csv")
transcript_stability <- transcript_stability[,c("transcript_id","gcsi_ccle_spearman","gcsi_gdsc_spearman","gdsc_ccle_spearman")]

gcsi_ccle <- transcript_stability[,c("transcript_id","gcsi_ccle_spearman")]
gcsi_gdsc <- transcript_stability[,c("transcript_id","gcsi_gdsc_spearman")]
gdsc_ccle <- transcript_stability[,c("transcript_id","gdsc_ccle_spearman")]

quantiles <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35)


# compute stable vs unstable overlap across ALL quantiles

for (q in quantiles){
  
  gcsi_ccle_stable <- gcsi_ccle[order(gcsi_ccle$gcsi_ccle_spearman, decreasing = TRUE),][1:floor(nrow(transcript_stability)*q),]
  gcsi_gdsc_stable <- gcsi_gdsc[order(gcsi_gdsc$gcsi_gdsc_spearman, decreasing = TRUE),][1:floor(nrow(transcript_stability)*q),]
  gdsc_ccle_stable <- gdsc_ccle[order(gdsc_ccle$gdsc_ccle_spearman, decreasing = TRUE),][1:floor(nrow(transcript_stability)*q),]
  
  gcsi_ccle_unstable <- gcsi_ccle[order(gcsi_ccle$gcsi_ccle_spearman, decreasing = FALSE),][1:floor(nrow(transcript_stability)*q),]
  gcsi_gdsc_unstable <- gcsi_gdsc[order(gcsi_gdsc$gcsi_gdsc_spearman, decreasing = FALSE),][1:floor(nrow(transcript_stability)*q),]
  gdsc_ccle_unstable <- gdsc_ccle[order(gdsc_ccle$gdsc_ccle_spearman, decreasing = FALSE),][1:floor(nrow(transcript_stability)*q),]
}
  
  
## Kuncheva Index (KI) - also known as consistency index ## --> STABILITY BY INDEX
quantiles <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35)

kuncheva_df <- data.frame(matrix(nrow=7, ncol = 2))
colnames(kuncheva_df) <- c("stable","unstable")
rownames(kuncheva_df) <- paste0("q",quantiles)

for (q in quantiles){
  
  gcsi_ccle_stable <- gcsi_ccle[order(gcsi_ccle$gcsi_ccle_spearman, decreasing = TRUE),][1:floor(nrow(transcript_stability)*q),]
  gcsi_gdsc_stable <- gcsi_gdsc[order(gcsi_gdsc$gcsi_gdsc_spearman, decreasing = TRUE),][1:floor(nrow(transcript_stability)*q),]
  gdsc_ccle_stable <- gdsc_ccle[order(gdsc_ccle$gdsc_ccle_spearman, decreasing = TRUE),][1:floor(nrow(transcript_stability)*q),]
  
  gcsi_ccle_unstable <- gcsi_ccle[order(gcsi_ccle$gcsi_ccle_spearman, decreasing = FALSE),][1:floor(nrow(transcript_stability)*q),]
  gcsi_gdsc_unstable <- gcsi_gdsc[order(gcsi_gdsc$gcsi_gdsc_spearman, decreasing = FALSE),][1:floor(nrow(transcript_stability)*q),]
  gdsc_ccle_unstable <- gdsc_ccle[order(gdsc_ccle$gdsc_ccle_spearman, decreasing = FALSE),][1:floor(nrow(transcript_stability)*q),]
  
  
  fsets <- list()
  fsets[["gcsi_ccle"]] <- gcsi_ccle_stable$transcript_id
  fsets[["gcsi_gdsc"]] <- gcsi_gdsc_stable$transcript_id
  fsets[["gdsc_ccle"]] <- gdsc_ccle_stable$transcript_id
  
  stable <- stab.fs(fsets=fsets, N=nrow(transcript_stability), method="kuncheva")
  
  fsets <- list()
  fsets[["gcsi_ccle"]] <- gcsi_ccle_unstable$transcript_id
  fsets[["gcsi_gdsc"]] <- gcsi_gdsc_unstable$transcript_id
  fsets[["gdsc_ccle"]] <- gdsc_ccle_unstable$transcript_id
  
  unstable <-  stab.fs(fsets=fsets, N=nrow(transcript_stability), method="kuncheva")
  
  
  kuncheva_df[ paste0("q",q), "stable"] <- stable
  kuncheva_df[ paste0("q",q), "unstable"] <- unstable
  
  
}
tt <- c("0.05","0.10","0.15","0.20","0.25","0.30","0.35")
kuncheva_df$quantile <-  tt
#kuncheva_df$quantile <-  gsub("q","",rownames(kuncheva_df))
kuncheva_m <- melt(kuncheva_df)

# plot line graph
png("../results/supplementary_figure9.png", width=150, height=100, units='mm', res = 600, pointsize=80)
ggplot(kuncheva_m, aes(x=quantile, y=value, group=variable), size = 3) +
  geom_line(aes(color=variable), size = 1)+
  geom_point(aes(color=variable)) + xlab("Quantile") + ylab("Kuncheva Index") + labs(color='Stability') + 
  scale_color_manual(values=c("#8B7B96", "#71A2B6")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),    # keep panel borders (removed from theme_classic())
        text = element_text(size = 17), 
        legend.key.size = unit(0.5, 'cm'),
        axis.text=element_text(size=15), 
        axis.text.x = element_text(size=15), 
        axis.text.y = element_text(size=15), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15))
dev.off()

