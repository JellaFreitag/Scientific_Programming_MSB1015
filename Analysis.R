

### Analysis

## PCA
# expression matrix: expr_only_filtered_probes
# row = samples, column = genes

expr_for_pca <- t(expr_only_filtered_probes) # transposing
pca_norm <- prcomp(expr_for_pca, center = TRUE, scale. = TRUE)

var_explained <- summary(pca_norm)$importance[2, 1:2] * 100
pc1_var <- round(var_explained[1], 1)
pc2_var <- round(var_explained[2], 1)

pc1_var
pc2_var
#summary(pca_norm)
#head(pca_norm$x


subtypes <- sub("SAMPLE \\d+\\s+", "", sample_names) 
subtypes <- sub("\\s+Replicate$", "", subtypes)
classes <- factor(subtypes)
classes2 <- factor(ifelse(classes == "PE", "PE", "NonPE"))


## DEG 
# expression matrix: expr_only_filtered_probes
# classes2 is already factorarized --> PE vs NonPE

design <- model.matrix(~0 + classes2)  # keine Intercept-Spalte
colnames(design) <- levels(classes2)
design

fit  <- lmFit(expr_only_filtered_probes, design)
cont <- makeContrasts(PE - NonPE, levels = design)
fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2)                     # Ergebnis in fit2 behalten
deg  <- topTable(fit2, adjust.method="BH", number=Inf)

# Signifikante DEGs (z.B. adj.P.Val < 0.1 und |logFC| > 1)
sig_deg <- deg[deg$adj.P.Val < 0.1 & abs(deg$logFC) > 1, ]


# how many genes up- or downregulated in PE compared to NonPE
table(sig_deg$logFC > 0)


# create matrices
deg_matrix <- as.matrix(deg) # all genes
sig_deg_matrix <- as.matrix(sig_deg) #  only significant genes

# create files
# setwd("")
write.csv(deg_matrix, "DEG_all_genes.csv")
write.csv(sig_deg_matrix, "DEG_significant_genes.csv")

