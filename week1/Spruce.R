# Load packages
## spruce package
library(spruce)
## Single cell code
library(Seurat)
# Bayesian modeling code
library(mvtnorm)
library(ggplot2)
library(BayesLogit)
library(QRM)
library(truncnorm)
library(patchwork)
library(SeuratData)
InstallData("stxBrain")
brain1 <- LoadData("stxBrain", type = "anterior1")
brain2 <- LoadData("stxBrain", type = "anterior2")
brain1 <- SCTransform(brain1, assay = "Spatial", verbose = FALSE)
brain2 <- SCTransform(brain2, assay = "Spatial", verbose = FALSE)
brain <- merge(brain1,brain2)
DefaultAssay(brain) <- "SCT"
VariableFeatures(brain) <- c(VariableFeatures(brain1),VariableFeatures(brain2))
brain <- RunPCA(brain)
PCA_8 <- brain@reductions$pca@cell.embeddings[,1:8]
# coordinate data
coords_df1 <- brain@images$anterior1@coordinates[,c("row","col")]
colnames(coords_df1) <- c("Y","X")

coords_df2 <- brain@images$anterior2@coordinates[,c("row","col")]
colnames(coords_df2) <- c("Y","X")
coords_df2$X <- coords_df2$X + 160

coords_df <- rbind(coords_df1,coords_df2)
w1 <- rep(1,nrow(PCA_8)) # intercept covariate
w2 <- as.numeric(brain$slice) - 1 # slice effect covariate
W <- cbind(w1,w2) # covariate matrix
fit_brain <- fit_mvn_PG_smooth(Y = PCA_8,
                               W = W,
                               coords_df = coords_df,
                               K = 6,
                               r = 3,
                               nsim = 2000,
                               burn = 1000)