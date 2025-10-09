#install.packages("readxl")
#install.packages("tidyverse")
#install.packages("BiocManager")
#BiocManager::install("limma")

library(readxl)
library(tidyverse)
library(limma)


# load data

setwd("C:/Users/jalli/Nextcloud/Systems Bio/Scientific Programming/EndoData")

Rdata <- read_excel("GSE141549_Non-normalized_data_GA_illumina_expression_platform_HumanHT-12.xls", col_names = FALSE, skip = 3)

ExprData <- Rdata %>% select(30:ncol(.)) # only expression columns + detection pval

colnames(ExprData) <- ExprData[1,]
view(table(colnames(ExprData)))
x <-  ExprData[1,]
sum(grepl("PP", x))
#---------------------------------------------------

# take only det p value columns:

sum(is.na(ExprData)) # check for NAs: none

# select only Detection Pval columns
Det_p_matrix <- ExprData %>% select(seq(2,ncol(.), by = 2)) %>% slice(-1) %>% sapply(as.numeric) %>% as.matrix()
class(Det_p_matrix)
dim(Det_p_matrix)

# find max value overall
max_p <- max(Det_p_matrix, na.rm = TRUE)
max_p
class(max_p)

# 10 largest numbers
sort(Det_p_matrix, decreasing = TRUE)[1:1000] # only 1e+05 --> too large?? outliers?

# vectorizing for visualization
det_p_vec <- as.numeric(as.vector(Det_p_matrix))
class(det_p_vec)


# bring detections scores to typical p-value ranges:
norm_Pvals <- Det_p_matrix/max_p

# how much significant? test for different thresholds
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
     xlab = "threshold", ylab = "significant data (%)", main = "Significant data vs. thresholds" )

good_probes <- norm_Pvals<0.05
table(good_probes)
1551760 /(1855496 + 1551760)*100

# only 45% are significant!? Yes! and it makes sense apparently for illumina microarray measurements!

# visualization

sigf_thresholds <- c(0.2, 0.15, 0.1, 0.08, 0.06, 0.05)

hist(norm_Pvals, breaks = 40, main = "Normalized p-values", xlab = "Scaled p-value")
abline(v = sigf_thresholds, col = "blue", lwd = 2, lty = 2)

#-----------------------------------------------------------------------------------------------------------------

######## Expression values ########
# extract expression data/remove p-values
expr_only <- ExprData %>% select(seq(1,ncol(.), by = 2)) %>% slice(-1) %>% sapply(as.numeric) %>% as.matrix()

#retreive GeneID vector
ILMN_Gene <- Rdata %>% select(6) %>% slice(-1)
ILMN_Gene <- ILMN_Gene$...6
ILMN_Gene

rownames(expr_only) <- ILMN_Gene

# visualization
raw <- boxplot(expr_only)
view(raw)

#Quality Control
sum(is.na(ExprData)) # check for NAs: none
sum(is.na(expr_only))

pca <- prcomp(t(expr_only), scale. = TRUE)
plot(pca$x[,1], pca$x[,2], main="PCA of raw data")
text(pca$x[,1], pca$x[,2], labels=colnames(expr_only), pos=1)

distMat <- dist(t(expr_only))
hc <- hclust(distMat)
plot(hc, main="Hierarchical clustering of samples")

##### --> probably no outliers

# densitiy plots?

# data distribution? can i see  the different classes?



######### p-values ##########


# filter out every gene that is insignificant
# expression data -> odd columns
# p-values -> even columns

#ExprData_filtered <- ExprData_scaled

#for (i in seq(1, ncol(ExprData_scaled), by = 2)) { # loop over all expression columns
#  expression_col <- i
#  p_column <- i+1
  
#  x <- ExprData_scaled[[p_column]] > 0.05 # check if p-value exceeds threshold, gives logic vector
#  ExprData_filtered[[expression_col]][x] <- NA
#  ExprData_filtered[[p_column]][x] <- NA
# }




######### p-values ##########
#### new approach: remove whole genes (rows), when average expression is below cutoff 

# scale p-values into range [0,1]

scale_Pvals <- function(ma){
  
  ma <- as.numeric(ma)
  ma / max(ma, na.rm = TRUE)
  
}


# select only Detection Pval columns
Det_p_matrix <- ExprData %>% select(seq(2,ncol(.), by = 2)) %>% slice(-1) %>% sapply(as.numeric) %>% as.matrix()
# scale p-value matrix
Det_p_matrix <- as.data.frame(Det_p_matrix)
p_vals_only <- Det_p_matrix %>% mutate(across(.cols = seq(1, ncol(.)), .fns = scale_Pvals))


### decide for threshold:
# logical matrix for each p-value; rowMeans counts all TRUE (1) and FALSE (0) per row together and gives back the mean
row_fraction <- rowMeans(p_vals_only < 0.1, na.rm = TRUE) # threshold = 0.1
row_fraction

### decide for cutoff:
# visual determination of best cutoff choice:
numb_genes <- length(row_fraction)
hist(row_fraction, breaks = 100, col = "orange")

# mean can be used to compare with cutoff: extract the 'good' genes
good_genes <- row_fraction >= 0.8
expr_only_filtered <- expr_only[good_genes,]
dim(expr_only_filtered) # 16818 genes left when 0.05 --> 26,22%
                        # 18568 genes left when 0.1 --> 28.18%
                  # Sagt ja nur aus, wie viele genes sich im average vom background noise abheben


## visualization
# non-filtered
plot(expr_only)
hist(expr_only)
# filtered
plot(expr_only_filtered)
hist(expr_only_filtered)


# log2-transformation
expr_only_log2 <- log2(expr_only_filtered+1)
plot(expr_only_log2)
hist(expr_only_log2)

# quantile normalization using limma
expr_only_norm <- normalizeBetweenArrays(expr_only_log2, method = "quantile")
plot(expr_only_norm)
hist(expr_only_norm)


# matrix with only gene expression values, vector with geneIDs
ILMN_Gene <- row.names(expr_only_norm)
table(ILMN_Gene)
expr_only_filtered_probes <- avereps(expr_only_norm, ID=ILMN_Gene)
          # --> leaves with 24.2% of the data
dim(expr_only_filtered_probes)


# transpose matrix
#DataExp <- as.data.frame(t(expr_only_filtered_probes))
#rownames(DataExp) <- paste0("Sample", 1:nrow(DataExp))
#dim(DataExp)
#class(DataExp)


# replicates --> how similar are they?
# retreive sample names + add to pre-processed matrix
sample_names <- as.character(ExprData[1,])
sample_names <- sample_names[seq(1, length(sample_names), by = 2)]

colnames(expr_only_filtered_probes) <- sample_names

# index of replicates + original
rep_idx <- grep("Replicate", colnames(expr_only_filtered_probes))
orig_idx <- rep_idx - 1 # original


# pearson and spearman correlation
cor_pearson <- sapply(seq_along(rep_idx), function(i) {
  cor(expr_only_filtered_probes[, orig_idx[i]], expr_only_filtered_probes[, rep_idx[i]], method = "pearson")
})
cor_pearson

cor_spearman <- sapply(seq_along(rep_idx), function(i) {
  cor(expr_only_filtered_probes[, orig_idx[i]], expr_only_filtered_probes[, rep_idx[i]], method = "spearman")
})
cor_spearman

# visualization
par(mfrow=c(1,2))  # 1 Reihe, 2 Spalten
hist(cor_pearson, main="Pearson Correlation", xlim=c(0,1), col="skyblue", breaks=10)
hist(cor_spearman, main="Spearman Correlation", xlim=c(0,1), col="salmon", breaks=10)
par(mfrow=c(1,1))





### bootstrapping
# data = expr_only_filtered_probes
dim(expr_only_filtered_probes)



# extract subtypes
library(stringr)
subtypes <- sub("SAMPLE \\d+\\s+", "", sample_names) 
subtypes <- sub("\\s+Replicate$", "", subtypes)
classes <- factor(subtypes)
table(classes)

# define classes: PE vs. the rest
classes2 <- factor(ifelse(classes == "PE", "PE", "NonPE"))
table(classes2) 

# actual bootstrapping
set.seed(99)

n_boot <- 1000
pe_idx <- which(classes2 == "PE")
nonpe_idx <- which(classes2 == "NonPE")

bootstrap_indices <- lapply(1:n_boot, function(i) {
  sample(nonpe_idx, length(nonpe_idx), replace = TRUE)  # 30 Non-PE, resampled
})
bootstrap_indices # indices of all 1000 possibilities of re-sampling, PE fix, NonPE flex


######################################### Analysis ##################################
 # create loop with analysis (PCA or DEG) to avoid memorizing all the 1000 different matrices


# test GitHub



