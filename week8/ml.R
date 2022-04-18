# BiocManager::install("Hiiragi2013")
library("Hiiragi2013")
library(dplyr)
library(GGally)
library(glmnet)
data("x")
probes = c("1426642_at", "1418765_at", "1418864_at", "1416564_at")
embryoCells = t(Biobase::exprs(x)[probes, ]) |> as_tibble() |>
  mutate(Embryonic.day = x$Embryonic.day) |>
  dplyr::filter(x$genotype == "WT")
annotation(x)
library("mouse4302.db")
anno = AnnotationDbi::select(mouse4302.db, keys = probes,
                             columns = c("SYMBOL", "GENENAME"))
anno
mt = match(anno$PROBEID, colnames(embryoCells))
colnames(embryoCells)[mt] = anno$SYMBOL
embryoCells
ggpairs(embryoCells, mapping = aes(col = Embryonic.day),
        columns = anno$SYMBOL, upper = list(continuous = "points"))

# lda
ec_lda = lda(Embryonic.day ~ Fn1 + Timd2 + Gata4 + Sox7,
             data = embryoCells)
round(ec_lda$scaling, 1)
ec_rot = predict(ec_lda)$x |> as_tibble() |>
  mutate(ed = embryoCells$Embryonic.day)
ec_lda2 = lda(ec_rot[, 1:2], predict(ec_lda)$class)

make1Dgrid = function(x) {
  rg = grDevices::extendrange(x)
  seq(from = rg[1], to = rg[2], length.out = 100)
}
ec_grid = with(ec_rot, expand.grid(
  LD1 = make1Dgrid(LD1),
  LD2 = make1Dgrid(LD2)))
ec_grid$edhat = predict(ec_lda2, newdata = ec_grid)$class
ggplot() +
  geom_point(aes(x = LD1, y = LD2, colour = ed), data = ec_rot) +
  geom_raster(aes(x = LD1, y = LD2, fill = edhat),
              data = ec_grid, alpha = 0.4, interpolate = TRUE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed()

# logistic regression
sx = x[, x$Embryonic.day == "E3.25"]
sx$genotype
embryoCellsClassifier = cv.glmnet(t(Biobase::exprs(sx)), sx$genotype,
                                  family = "binomial", type.measure = "class")
plot(embryoCellsClassifier)

# SVM
library(e1071)
svm = svm(t(Biobase::exprs(sx)), sx$genotype)
pred <- predict(svm, t(Biobase::exprs(sx)))
table(sx$genotype, pred)
