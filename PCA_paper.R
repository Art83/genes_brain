library(PerformanceAnalytics)
library(factoextra)
library(FactoMineR)
library(corrplot)
library(xlsx)
library(dplyr)
library(rgl)



dat_30 <- read.xlsx("data/Summary_BLA_genes_30.xlsx", sheetIndex = 4, header = T)
dat_90 <- read.xlsx("data/PCR_summary_BLA.xlsx", sheetIndex = 18, as.data.frame = T, header = T)
dat_270 <- read.xlsx("data/Summary_BLA_genes_30.xlsx", sheetIndex = 6, header = T)


cond1 <- c("second-order", 'first-order')
dat_30_goi <- filter(dat_30, subgroup %in% cond1)
dat_90_goi <- filter(dat_90, subgroup %in% cond1)
dat_30_goi <- dat_30_goi[-11, ]
dat_90_goi <- dat_90_goi[-c(),]
mydata <- data.frame(scale(dat_30_goi[,-c(1,2)]))
mydata <- data.frame(scale(dat_90_goi[,-c(1,2)]))
mydata <- data.frame(dat_90_goi[,-c(1,2)])
mydata <- data.frame(dat_30_goi[,-c(1,2)])
res.pca$ind$coord



res.pca <- PCA(mydata, graph = F)


pc <- prcomp(mydata)
comp <- data.frame(pc$x[,1:4])
plot(comp, pch=16, col=rgb(0,0,0,0.5))
plot3d(comp[,1], comp[,2], comp[,3])
k <- kmeans(comp, 2, nstart=25, iter.max=1000)

clust <- c(rep(1, 8), rep(2, 6))

library(RColorBrewer)
library(scales)
palette(alpha(brewer.pal(9,'Set1'), 0.5))
plot(comp, col=clust, pch=16)

plot3d(comp$PC1, comp$PC2, comp$PC3, col=clust, size = 14)
plot3d(comp$PC1, comp$PC2, comp$PC4, col=clust, size = 14)


eig.values <- get_eigenvalue(res.pca)

fviz_eig(res.pca, addlabels = T, ylim = c(0,90)) # Scree plot , PC1:4??



open3d()
plot3d(comp$PC1, comp$PC2, comp$PC3, col=clust, size = 14)

rgl.postscript("persp3dd.pdf","pdf")

fviz_pca_var(res.pca, col.var = 'black')

fviz_contrib(res.pca, choice = 'var', axes = 1:3, top = 10)




mydata <- scale(participants_table[-1, seq(10, 24, 2)])
rownames(mydata) <- participants_table$PIN[-1]
wss <- (nrow(mydata) - 1) * sum(apply(mydata, 2, var))

for (i in 2:15){
  wss[i] <- sum(kmeans(mydata, centers = i)$withinss)
}

plot(1:15, wss, type='b', xlab = 'Number of clusters', ylab='Within groups sum of squares')


km.out <- kmeans(mydata, 2)

km.out$cluster


plot(mydata, col=(km.out$cluster+1), main = 'Clustering', xlab='', ylab='', pch=20, cex=2)



participants_table$PIN

fviz_nbclust(mydata, kmeans)



fviz_cluster(km.out, data = mydata)+
  theme_classic()


fviz_pca_biplot(res.pca,
                geom.ind = 'point',
                pointshape = 21,
                pointsize = 2.5,
                col.ind = 'blue',
                repel = T)

fviz_pca_biplot(res.pca, repel = T,
                col.var = 'red',
                col.ind = 'blue',
                axes = c(1,2))


fviz_pca_biplot(res.pca, 
                # Fill individuals by groups
                geom.ind = "point",
                pointshape = 21,
                pointsize = 2.5,
                fill.ind = dat_30_goi$subgroup,
                #col.ind = dat_90_goi$subgroup,
                # Color variable by groups
                legend.title = list(fill = "Species"),
                repel = TRUE,        # Avoid label overplotting
                axes = c(1,2)
)+
  ggpubr::fill_palette("jco")+      # Indiviual fill color
  ggpubr::color_palette("npg")      # Variable colors


iris.pca <- PCA(iris[,-5], graph = FALSE)
iris$Species
fviz_pca_ind(iris.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = iris$Species, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)

dat_90_goi$subgroup
