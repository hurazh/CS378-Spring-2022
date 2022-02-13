BiocManager::install("GEOquery")
BiocManager::install("preprocessCore")
install.packages('pheatmap')
install.packages('matrixStats')

library(GEOquery)
library(preprocessCore)
library(matrixStats)
library(pheatmap)

file <- list.files(pattern="GSE2125")
if (length(file) == 0) {
  geo <- getGEO("GSE2125", destdir=".")
  e <- geo[[1]]
} else {
  e <- getGEO(filename=file)
}

head(exprs(e))
head(pData(e))
head(fData(e))
e$extract_protocol_ch1[1]
e$characteristics_ch1
table(e$characteristics_ch1)

range(exprs(e))
boxplot(exprs(e),range=0)

exprs(e) <- normalize.quantiles(exprs(e))
boxplot(exprs(e),range=0)

exp_val <- rowVars(exprs(e))
o <- head(order(exp_val, decreasing=TRUE),200)
pheatmap(exprs(e)[o,],
         annotation_col=pData(e)["characteristics_ch1"],
         show_rownames=FALSE, show_colnames=FALSE)

