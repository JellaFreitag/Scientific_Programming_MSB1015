## visualizations


############ add QC plots, threshold plots etc. !!





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

# check expression distribution again after processing
boxplot(expr_only_filtered_probes, outline = FALSE, las = 2, main = "Expression distribution per sample",
        ylab = "log2(Expression)", xaxt = "n", xlab = "samples", col = "lightblue")



# sample relationship
par(mfrow=c(1,2),oma = c(0, 0, 3, 0))  
hist(cor_pearson, main="Pearson Correlation", xlim=c(0,1), xlab = "Correlation coefficient", ylab = "Number of replicate pairs", 
     col="skyblue", breaks=10)
hist(cor_spearman, main="Spearman Correlation", xlim=c(0,1), xlab = "Correlation coefficient", ylab = "", 
     col="yellow", breaks=10)
title("Distribution of correlation coefficients between replicates", outer = TRUE)
par(mfrow=c(1,1))
par(oma = c(0, 0, 0, 0))


# class imbalance
subtypes <- sub("SAMPLE \\d+\\s+", "", sample_names) 
subtypes <- sub("\\s+Replicate$", "", subtypes)
classes <- factor(subtypes)
table(classes)

# all classes
barplot(table(classes), col = "purple4", main = "Sample distribution across subtypes", xlab = "Classes", ylab = "Number of samples")

# PE vs. the rest
classes2 <- factor(ifelse(classes == "PE", "PE", "NonPE"))
table(classes2) 

barplot(table(classes2), col = "purple4", main = "Sample distribution across two classes", xlab = "Classes", ylab = "Number of samples")


## Analysis

# PCA

# PE vs. NonPE
classes2
colors <- c("PE" = "blue", "NonPE" = "green")
point_colors <- colors[classes2]

plot(pca_norm$x[,1], pca_norm$x[,2],
     xlab = paste0("PC1 (", pc1_var, "% variance)"), ylab =  paste0("PC2 (", pc2_var, "% variance)"),
     main = "PCA of gene expression in endometriosis samples", col = point_colors, pch = 20)
legend("topright", legend = names(colors), col = colors, pch = 20, x.intersp = 0.6,
       y.intersp = 0.8, bty = "n", bg = "transparent")


# all subtypes
classes
group_colors <- c("DiEIn" = "red",
                  "PE"    = "blue",
                  "PeLB"  = "green",
                  "PeLR"  = "orange",
                  "PeLW"  = "purple",
                  "PP"    = "brown",
                  "SuL"   = "pink")

point_colors2 <- group_colors[classes]


plot(pca_norm$x[,1], pca_norm$x[,2],
     xlab = paste0("PC1 (", pc1_var, "% variance)"), ylab =  paste0("PC2 (", pc2_var, "% variance)"),
     main = "PCA of gene expression in endometriosis subtypes",
     col = point_colors2, pch = 20)

legend("topright", legend = names(group_colors), col = group_colors, pch = 20, x.intersp = 0.6,
       y.intersp = 0.8, bty = "n", bg = "transparent")
#text(pca_norm$x[,1], pca_norm$x[,2], labels = rownames(expr_for_pca), pos = 3, cex = 0.7)



# visualize DEG: Volcano Plot
plot(deg$logFC, -log10(deg$adj.P.Val),
     pch = 19, cex = 0.5,
     col = ifelse(deg$adj.P.Val < 0.1 & abs(deg$logFC) > 1, "blue", "grey"),
     xlab = "log2 Fold Change", ylab = "-log10 adjusted p-value",
     main = "Volcano plot of DEGs (PE vs NonPE)")




