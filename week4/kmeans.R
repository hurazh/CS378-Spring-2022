setwd("~/Desktop/compbio/week4")
load(file='tissuesGeneExpression.rda')
exp <- data.frame(t(e))

# KMeans
set.seed(1)
# first 2 gene
km_1 <- kmeans(exp[,1:2], centers=7)
names(km_1)
ggplot(exp, aes(X1007_s_at, X1053_at, color=tissue)) +   
  geom_point(size=2, pch=16)
cluster_1 <- as.character(km_1$cluster)
ggplot(exp, aes(X1007_s_at, X1053_at, color=cluster_1)) +   
  geom_point(size=2, pch=16)
table(true=tissue,cluster=km_1$cluster)
# all genes
km_2 <- kmeans(exp, centers=7)
table(true=tissue,cluster=km_2$cluster)
d <- dist(exp)
mds <- cmdscale(d)
distance <- data.frame(X=mds[,1],
                       Y=mds[,2])
ggplot(distance, aes(X, Y, color=tissue)) +   
       geom_point(size=2, pch=16)
cluster_2 <- as.character(km_2$cluster)
ggplot(distance, aes(X, Y, color=cluster_2)) +   
       geom_point(size=2, pch=16)

# hierarchical clustering
hc <- hclust(d)
hc
plot(hc, labels=tissue, cex=0.5)
hcluster_h <- cutree(hc, h=120)
table(true=tissue, cluster=hcluster_h)
hcluster_k <- cutree(hc, k=7)
table(true=tissue, cluster=hcluster_k)
