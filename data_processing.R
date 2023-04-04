
if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (! requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery")
}
if (!requireNamespace("biomaRt", quietly = TRUE, force = TRUE)) {
  BiocManager::install("biomaRt")
}
library(GEOquery)
library(knitr)
library(edgeR)
library(biomaRt)

# downloading GEO dataset
GSE <- getGEO("GSE221786", GSEMatrix=FALSE)

# Information about Platform
GPL <- names(GPLList(GSE))[1]
GPL_info <- Meta(getGEO(GPL))

# get the expression data
sfiles <- getGEOSuppFiles('GSE221786')
fnames <- rownames(sfiles)
count_exp = read.delim(fnames[1],header=TRUE, check.names = FALSE)
samples <- data.frame(lapply(colnames(count_exp)[3:30],
                             FUN=function(x){unlist(strsplit(x,
                                                             split = "\\_"))[c(1,2)]}))
colnames(samples) <- colnames(count_exp)[3:30]
rownames(samples) <- c("group","gender")
samples <- data.frame(t(samples))

cpms <- edgeR::cpm(count_exp[,3:30])
rownames(cpms) <- count_exp[,1]
# Get rid of low counts, there are 2 samples of each group so n=2
keep = rowSums(cpms > 1) >= 2
count_exp_filtered <- count_exp[keep, ]
# After filtering out genes that have low counts, there are 15037 genes remaining
nrow(count_exp) - nrow(count_exp_filtered) # 10191

# Get the summarized counts for each gene
summarized_gene_counts_filtered <- sort(table(count_exp_filtered$gene),decreasing = TRUE)
length(summarized_gene_counts_filtered[which(summarized_gene_counts_filtered>1)])
# The problem of duplicated genes has been solved partially
kable(summarized_gene_counts_filtered[
  which(summarized_gene_counts_filtered>1)[1:10]],
  format="html")
count_exp_filtered <- count_exp_filtered[!duplicated(count_exp_filtered$gene),]
rownames(count_exp_filtered) <- count_exp_filtered[,1]

filtered_data_matrix <- as.matrix(count_exp_filtered[,3:30])
# Rownames be the ensembl id
rownames(filtered_data_matrix) <- count_exp_filtered$ensembl_id
d = DGEList(counts=filtered_data_matrix, group=samples$group)
# Normalization
d = calcNormFactors(d)
normalized_counts <- cpm(d)
rownames(normalized_counts) <- count_exp_filtered$gene
nrow(normalized_counts)
# HUGO symbols as rownames of the dataframe
kable(normalized_counts[1:5,1:5], format = "html")
# normalized_counts is our result data frame. Use for A2.
