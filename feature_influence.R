suppressMessages(library(caret))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
set.seed(200)  

#Read in transcript stability data generated from ../code/transcript_stability.R
transcript_stability <- read.csv("../results/transcript_stability.csv")


#Make a copy for each dataset pair
transcript_stability_gcsi_ccle <- transcript_stability
transcript_stability_gdsc_ccle <- transcript_stability
transcript_stability_gcsi_gdsc <- transcript_stability


#Isolate MOST stable/unstable transcripts (15% quantile)
#gCSI/CCLE
stable <- transcript_stability_gcsi_ccle[order(-transcript_stability_gcsi_ccle$gcsi_ccle_spearman),]
stable <- stable[1:9037,] #top 15%
unstable <- transcript_stability_gcsi_ccle[order(transcript_stability_gcsi_ccle$gcsi_ccle_spearman),]
unstable <- unstable[1:9037,]
transcript_stability_gcsi_ccle <- rbind(stable,unstable)

#gCSI/GDSC
stable <- transcript_stability_gcsi_gdsc[order(-transcript_stability_gcsi_gdsc$gcsi_gdsc_spearman),]
stable <- stable[1:9037,]
unstable <- transcript_stability_gcsi_gdsc[order(transcript_stability_gcsi_gdsc$gcsi_gdsc_spearman),]
unstable <- unstable[1:9037,]
transcript_stability_gcsi_gdsc <- rbind(stable,unstable)

#GDSC/CCLE
stable <- transcript_stability_gdsc_ccle[order(-transcript_stability_gdsc_ccle$gdsc_ccle_spearman),]
stable <- stable[1:9037,]
unstable <- transcript_stability_gdsc_ccle[order(transcript_stability_gdsc_ccle$gdsc_ccle_spearman),]
unstable <- unstable[1:9037,]
transcript_stability_gdsc_ccle <- rbind(stable,unstable)


#Extract features + stability from df
transcript_stability_gcsi_ccle <- transcript_stability_gcsi_ccle[,c('gdsc_median','gc','n_exon','single_read_mappability','multi_read_mappability','length','three_prime_bias','five_prime_bias','gcsi_ccle_spearman')]
transcript_stability_gcsi_gdsc <- transcript_stability_gcsi_gdsc[,c('ccle_median','gc','n_exon','single_read_mappability','multi_read_mappability','length','three_prime_bias','five_prime_bias','gcsi_gdsc_spearman')]
transcript_stability_gdsc_ccle <- transcript_stability_gdsc_ccle[,c('gcsi_median','gc','n_exon','single_read_mappability','multi_read_mappability','length','three_prime_bias','five_prime_bias','gdsc_ccle_spearman')]

#shuffle data
transcript_stability_gcsi_ccle <- transcript_stability_gcsi_ccle[sample(1:nrow(transcript_stability_gcsi_ccle)),]
transcript_stability_gcsi_gdsc <- transcript_stability_gcsi_gdsc[sample(1:nrow(transcript_stability_gcsi_gdsc)),]
transcript_stability_gdsc_ccle <- transcript_stability_gdsc_ccle[sample(1:nrow(transcript_stability_gdsc_ccle)),]


########################################################
########### Permutation Feature Importance #############
########################################################


################# 
# Reduce size of data.frame for Code Ocean, comment out following 3 lines for full result
transcript_stability_gcsi_ccle <- transcript_stability_gcsi_ccle[1:500,]
transcript_stability_gcsi_gdsc <- transcript_stability_gcsi_gdsc[1:500,]
transcript_stability_gdsc_ccle <- transcript_stability_gdsc_ccle[1:500,]
#################


## Run Linear Regression
runLinearReg <- function(x, label) {
  
  # Create dataframe to store feature scores
  feature_score <- data.frame(matrix(ncol=4, nrow = 0))
  colnames(feature_score) <- c("feature", "baseline", "permuted", "dataset")
  
  # Change colnames of inputted dataframe
  colnames(x) <- c("median", colnames(x)[2:8], "spearman")
  
  #10-fold CV
  fitControl <- trainControl(method = "repeatedcv",
                             number = 10,
                             repeats = 5) #change to 10 for full results
  
  
  # Loop through each feature
  for (i in 1:length(colnames(x)[1:8])) {
    
    feature <- colnames(x)[i]
    
    # get baseline score
    gbmFit <- train(spearman ~ ., data = x, 
                    method = "lm", 
                    trControl = fitControl,
                    verbose = FALSE)
    #Save MSE score
    baseline_score <- (gbmFit$results$RMSE)^2
    
    # Permute each feature 20,000 times
    for (k in 1:2) { # change to 20,000 for full results
      
      # shuffle feature column
      shuffle_feat <- sample(x[,i])
      x[,i] <- shuffle_feat
      
      gbmFit <- train(spearman ~ ., data = x, 
                      method = "lm", 
                      trControl = fitControl,
                      verbose = FALSE)
      
      # Record permuted MSE score
      permutation_score <- (gbmFit$results$RMSE)^2
      
      score <- data.frame(feature = feature, baseline = baseline_score, permuted = permutation_score, dataset = label)
      feature_score <- rbind(feature_score, score)
    }
  }
  return(feature_score)
}


lm_gcsi_ccle <- runLinearReg(transcript_stability_gcsi_ccle, "gcsi_ccle")
lm_gcsi_gdsc <- runLinearReg(transcript_stability_gcsi_gdsc, "gcsi_gdsc")
lm_gdsc_ccle <- runLinearReg(transcript_stability_gdsc_ccle, "gdsc_ccle")

lm_feature_results <- rbind(lm_gcsi_ccle, lm_gcsi_gdsc, lm_gdsc_ccle)



###############
#### PLOTS ####
###############

#rename datasets
lm_feature_results <- lm_feature_results %>% 
  mutate(dataset = recode(dataset, "gcsi_ccle" = "gCSI/CCLE", "gcsi_gdsc" = "gCSI/GDSC", "gdsc_ccle" = "GDSC/CCLE"))

#rename features
lm_feature_results <- lm_feature_results %>% 
  mutate(feature = recode(feature, "median" = "Median Expression", "gc" = "GC%", "n_exon" = "Number of Exons", "length" = "Length",
                          "single_read_mappability" = "Single Read Mappability", "multi_read_mappability" = "Multi Read Mappability",
                          "three_prime_bias" = "Three Prime Bias", "five_prime_bias" = "Five Prime Bias"))


lm_feature_results$factor <- factor(lm_feature_results$feature, 
                                    levels=c('Median Expression', 'GC%', 'Number of Exons', 'Length', 'Single Read Mappability', 'Multi Read Mappability', 'Three Prime Bias', 'Five Prime Bias'))

### Model Analysis ###

lm_model <- melt(lm_feature_results)


p1 <- ggplot(data = lm_model) +
  geom_density(aes(x = value, fill = variable), color = "black", alpha = 0.5) +
  facet_wrap(factor ~ ., ncol = 4) +
  scale_fill_manual(values=c("#71A2B6", "#982649")) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),    
        text = element_text(size = 15), 
        legend.key.size = unit(1, 'cm'),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 17), 
        axis.text.x = element_text(size=12, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Linear Regression Permutation Results", 
       y = "Density", x = "MSE") +
  guides(fill = guide_legend(title = "Score", ncol = 2))

png("../results/supplementary_figure10.png", width=275, height=175, units='mm', res = 600, pointsize=80)
p1
dev.off()




#### Feature Influence ####

lm_feature_results$difference <- with(lm_feature_results, baseline - permuted)


p1 <- ggplot(data = lm_feature_results, aes(x=feature, y=difference)) + geom_boxplot(aes(color=dataset), lwd=0.75) +  
  ylab(expression('MSE - MSE'['permuted'])) + xlab("\nFeature") + 
  scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("gCSI/CCLE","gCSI/GDSC","GDSC/CCLE"), values = c("#392C57", "#51C7AD", "#3670A0")) + 
  scale_x_discrete(limits = c('Median Expression', 'GC%', 'Number of Exons', 'Length', 'Single Read Mappability', 'Multi Read Mappability', 'Three Prime Bias', 'Five Prime Bias'),
                   labels=c('Median\nExpression', 'GC%', 'Number\nof Exons', 'Length', 'Single Read\nMappability','Multi Read\nMappability', 'Three Prime\nBias','Five Prime\nBias')) + 
  geom_hline(yintercept = 0, size = 0.5, linetype = "dotted", color = "black") +
  theme_classic() + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        text = element_text(size = 15), 
        legend.key.size = unit(0.6, 'cm'),
        axis.text.x = element_text(size=12, vjust = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 17), 
        axis.text.y = element_text(size=15), 
        legend.text = element_text(size=15), 
        legend.position = "bottom") + 
  labs(title = "Linear Regression Permutation Feature Influence") +
  guides(guide_legend(title.vjust = 2.8)) 


png("../results/figure3.png", width=275, height=150, units='mm', res = 600, pointsize=80)
p1
dev.off()