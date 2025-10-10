# make sure that the variables and matrices from Script1 are loaded!

############### Data filtering for significance ############### 

### first approach: filter out every element from ExprData matrix
# every expression value and the corresponding p-value

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

### problem here: there is a lot of missing values in the matrix, which hinders further processing and analysis



### new approach: remove whole genes (rows), when fraction expression is below cutoff 
# expression matrix: ExprData

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


## decide for threshold:
# logical matrix for each p-value; rowMeans counts all TRUE (1) and FALSE (0) per row together 
# gives back the mean of all the elements for each gene, for p-values below 0.1
row_fraction <- rowMeans(p_vals_only < 0.1, na.rm = TRUE) # threshold = 0.1
row_fraction

## decide for fraction cutoff:
# visual determination of best cutoff choice:
numb_genes <- length(row_fraction)
hist(row_fraction, breaks = 100, col = "orange", main = "Distribution of gene significance using row fraction (p < 0.1)", 
     ylab = "Number of genes", xlab = "Row fraction")

# mean can be used to compare with cutoff: extract the 'good' genes
good_genes <- row_fraction >= 0.8
expr_only_filtered <- expr_only[good_genes,]
dim(expr_only_filtered) # 16818 genes left when 0.05 --> 26,22%
                        # 18568 genes left when 0.1 --> 28.18%
                        # Sagt ja nur aus, wie viele genes sich im average vom background noise abheben

# log2-transformation
expr_only_log2 <- log2(expr_only_filtered+1)

# quantile normalization using limma
expr_only_norm <- normalizeBetweenArrays(expr_only_log2, method = "quantile")



# filter for probes
# matrix with only gene expression values, vector with geneIDs
ILMN_Gene <- row.names(expr_only_norm)
table(ILMN_Gene)
expr_only_filtered_probes <- avereps(expr_only_norm, ID=ILMN_Gene) # average replicates
# --> leaves with 24.2% of the data
dim(expr_only_filtered_probes)



# replicates of samples
# replicates --> how similar are the samples?
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
par(mfrow=c(1,2),oma = c(0, 0, 3, 0))  
hist(cor_pearson, main="Pearson Correlation", xlim=c(0,1), xlab = "Correlation coefficient", ylab = "Number of replicate pairs", 
     col="skyblue", breaks=10)
hist(cor_spearman, main="Spearman Correlation", xlim=c(0,1), xlab = "Correlation coefficient", ylab = "", 
     col="yellow", breaks=10)
title("Distribution of correlation coefficients between replicates", outer = TRUE)
par(mfrow=c(1,1))





