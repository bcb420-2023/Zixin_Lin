# A script that performs the QLF testing and gene ranking in A2
# but without R Markdown documentation. Intended to be sourced in A2_LukeZhang.Rmd.

# The result of this script is a sorted list of differentially expressed genes, with
# upregulated genes at the top and downregulated genes at the bottom.

# Download packages ===========================================================

if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (! requireNamespace("edgeR", quietly = TRUE)) {
  BiocManager::install("edgeR")
}
library(edgeR)


# construct model
model_design_pat <- model.matrix(~ samples$group+ samples$gender)

d = DGEList(counts=normalized_counts, group=samples$gender)
# estimate dispersion
d <- estimateDisp(d, model_design_pat)
# calculate normalization factors
d <- calcNormFactors(d)
# fit model to data
fit <- glmQLFit(d, model_design_pat)
# calculate differential expression
qlf <- glmQLFTest(fit)

# Get all the results
qlf_output_hits <- topTags(qlf, 
                           sort.by = "PValue", 
                           n = nrow(normalized_counts))
ranked_genes <- qlf_output_hits

