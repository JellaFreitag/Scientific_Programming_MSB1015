#Load packages
#install.packages("readxl")
#install.packages("tidyverse")
#install.packages("BiocManager")
#BiocManager::install("limma")

library(readxl)
library(tidyverse)
library(limma)
library(stringr)


# Load data
setwd("C:/Users/jalli/Nextcloud/Systems Bio/Scientific Programming/EndoData")

Rdata <- read_excel("GSE141549_Non-normalized_data_GA_illumina_expression_platform_HumanHT-12.xls", col_names = FALSE, skip = 3)

ExprData <- Rdata %>% select(30:ncol(.)) # only expression columns + detection pval



# Exploring data 
######### detection p-vals #########

# select only Detection Pval columns
Det_p_matrix <- ExprData %>% select(seq(2,ncol(.), by = 2)) %>% slice(-1) %>% sapply(as.numeric) %>% as.matrix()
class(Det_p_matrix)
dim(Det_p_matrix)


# how to deal with detection scores? --> scaling to common range --> what is max value?
# 1000 largest numbers
sort(Det_p_matrix, decreasing = TRUE)[1:1000] # only 1e+05 --> too large?? outliers?

# find max value overall
max_p <- max(Det_p_matrix, na.rm = TRUE)
max_p
class(max_p)

# vectorizing for visualization: densitiy plot
det_p_vec <- as.numeric(as.vector(Det_p_matrix))
class(det_p_vec)

dens <- density(det_p_vec)
plot(dens$x, dens$y, type = "l", main = "Density of detection scores", xlab = "Detection score", ylab = "Density")
#hist(det_p_vec)

# bring detections scores to typical p-value ranges:
norm_Pvals <- Det_p_matrix/max_p


# testing for different p-value threshold
sigf_thresholds <- c(0.2, 0.15, 0.1, 0.08, 0.06, 0.05, 0.01, 0.001)
sigf_data <- numeric(length(sigf_thresholds))

for (i in seq(sigf_thresholds)) {
  good_probes <- norm_Pvals<sigf_thresholds[i]
  num_true <- sum(good_probes)
  num_false <- sum(!good_probes)
  sigf_data[i] <- num_true/(num_false+num_true)*100
  
}
sigf_data


plot(sigf_thresholds,sigf_data, type = "b", ylim = c(0,100),pch=20, col = "blue", 
     xlab = "Threshold", ylab = "Significant data (%)", main = "Threshold effect on significance")


hist(norm_Pvals, breaks = 40, main = "Distribution of scaled detection scores with thresholds", xlab = "Scaled p-values")
abline(v = sigf_thresholds, col = "blue", lwd = 2, lty = 2)





######## Expression values ########
# extract expression data/remove p-values
expr_only <- ExprData %>% select(seq(1,ncol(.), by = 2)) %>% slice(-1) %>% sapply(as.numeric) %>% as.matrix()

#retreive GeneID vector
ILMN_Gene <- Rdata %>% select(6) %>% slice(-1)
ILMN_Gene <- ILMN_Gene$...6
ILMN_Gene

rownames(expr_only) <- ILMN_Gene

# visualization raw data
#raw <- boxplot(expr_only) # cannot really see much
#view(raw)


# Quality Control

sum(is.na(ExprData)) # check for NAs: none
sum(is.na(expr_only))

# PCA raw data - all classes
# retreive sample names 
sample_names <- as.character(ExprData[1,])
sample_names <- sample_names[seq(1, length(sample_names), by = 2)]
subtypes <- sub("SAMPLE \\d+\\s+", "", sample_names) 
subtypes <- sub("\\s+Replicate$", "", subtypes)
classes <- factor(subtypes)
classes
group_colors <- c("DiEIn" = "red",
                  "PE"    = "blue",
                  "PeLB"  = "green",
                  "PeLR"  = "orange",
                  "PeLW"  = "purple",
                  "PP"    = "brown",
                  "SuL"   = "pink")

all_class_colors <- group_colors[classes]

pca <- prcomp(t(expr_only), scale. = TRUE)
plot(pca$x[,1], pca$x[,2], xlab = "PC1", ylab = "PC2", main="PCA of raw data", col = all_class_colors, pch = 20)
legend("topright", legend = names(group_colors), col = group_colors, pch = 20, x.intersp = 0.6,
       y.intersp = 0.6, bty = "n", bg = "transparent")
text(pca$x[,1], pca$x[,2], labels=colnames(expr_only), pos=2)


# Hierarchical clustering
distMat <- dist(t(expr_only))
hc <- hclust(distMat)
plot(hc, main="Hierarchical clustering of samples - raw data", sub = "", xlab = "")

## --> color-code dendrogram would be nice..
## --> but no ouliers

boxplot(log2(expr_only + 1), outline = FALSE, las = 2, main = "Expression distribution per sample",
        ylab = "log2(Expression)", xaxt = "n", xlab = "samples", col = "lightblue")

