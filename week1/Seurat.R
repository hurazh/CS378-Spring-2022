brain1 <- LoadData("stxBrain", type = "anterior1")
brain1 <- SCTransform(brain1, assay = "Spatial", verbose = FALSE)
DefaultAssay(brain1) <- "SCT"
brain1@assays
brain1 <- RunPCA(brain1)
PCA_8 <- brain1@reductions$pca@cell.embeddings[,1:8]


zebrafish.seurat <- as.Seurat(zebrafish, counts = "counts",data = "logcounts")
zebrafish.seurat <- SCTransform(zebrafish.seurat, assay = "originalexp", verbose = FALSE)
DefaultAssay(zebrafish.seurat) <- "SCT"
zebrafish.seurat <- RunPCA(zebrafish.seurat)
coords_zebrafish <- colData_zebrafish
colnames(coords_zebrafish) <- c("Y","X")
PCA_8 <- zebrafish.seurat@reductions$pca@cell.embeddings[,1:8]
w1 <- rep(1,nrow(PCA_8))
w2 <- rep(0,nrow(PCA_8))
PCA_8 = PCA_8 + runif(dim(PCA_8)) * 1e-4
W <- cbind(w1,w2) # covariate matrix
fit_semisyn_27 <- fit_mvn_PG_smooth(Y = PCA_8,
                                   W = W,
                                   coords_df = coords_zebrafish,
                                   K = 6,
                                   nsim = 2000,
                                   burn = 1000)
get_scores(fit_semisyn_27)
fit_zebrafish_7 <- fit_mvn_PG_smooth(Y = PCA_8,
                                    W = W,
                                    coords_df = coords_zebrafish,
                                    K = 7,
                                    nsim = 2000,
                                    burn = 1000)
get_scores(fit_zebrafish_7)
fit_zebrafish_8 <- fit_mvn_PG_smooth(Y = PCA_8,
                                     W = W,
                                     coords_df = coords_zebrafish,
                                     K = 8,
                                     nsim = 2000,
                                     burn = 1000)
get_scores(fit_zebrafish_8)
coords_df <- data.frame(X = coords_zebrafish$X, Y = coords_zebrafish$Y)
hex_plt <- plot_spatial_hex(coords_zebrafish, as.integer(fit_zebrafish$z))
plot(hex_plt)
write.csv(fit_zebrafish_8$z, 'spruce_8.cvs')

