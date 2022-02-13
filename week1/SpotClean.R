library(SpotClean)

zebrafish_sim_mat <- Read10xRaw(count_dir = "/Users/hrzhang/Documents/projects/SpatialTranscriptomics/SGE/python/bleeding_sim/anisotropic/raw_feature_bc_matrix")
object.size(zebrafish_sim_mat)
zebrafish_info <- Read10xSlide(tissue_csv_file="/Users/hrzhang/Documents/projects/SpatialTranscriptomics/SGE/python/bleeding_sim/anisotropic/spatial/tissue_positions_list.csv")
sim_obj <- CreateSlide(count_mat = zebrafish_sim_mat, 
                          slide_info = zebrafish_info)
decont_aniso_obj <- SpotClean(sim_obj)
full_mat = as.matrix(decont_aniso_obj@assays@data@listData$decont)
write.csv(full_mat, 'anisotropic_spotclean_2_4.csv')

decont_sim_obj <- SpotClean(sim_obj)
full_mat = as.matrix(decont_sim_obj@assays@data@listData$decont)
write.csv(full_mat, 'gaussian_spotclean_2_3.csv')
row_sim = decont_sim_obj@metadata$slide$row
col_sim = decont_sim_obj@metadata$slide$col
write.csv(row_sim, 'row_sim.csv')
write.csv(col_sim, 'col_sim.csv')

which(sim_obj@assays@data@listData$raw@Dimnames[[1]] == "BRAFhuman")
