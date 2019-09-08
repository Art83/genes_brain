
library(xlsx)
library(PerformanceAnalytics)
library(factoextra)
library(ggbiplot)
library(ggplot2)
library(dplyr)

setwd("/home/arthur/Documents/R_learning")

# reading data
dat_30 <- read.xlsx("Summary_BLA_genes_30.xlsx", sheetIndex = 4, header = T)
dat_90 <- read.xlsx("PCR_summary_BLA.xlsx", sheetIndex = 18, as.data.frame = T, header = T)
dat_270 <- read.xlsx("Summary_BLA_genes_30.xlsx", sheetIndex = 6, header = T)

dat <- rbind(dat_30, dat_90, dat_270)

# handling 30 minutes dataset
cond2 <- c("second-order", 'first-order', "spc")
cond1 <- c("second-order", 'first-order')
cond3 <- c('first-order', "spc")
cond4 <- c("second-order", "spc")


dat_30_goi <- filter(dat_30, subgroup %in% cond1)
dat_30_goi <- filter(dat_30, subgroup %in% cond2)
dat_30_goi <- filter(dat_30, subgroup %in% cond3)
dat_30_goi <- filter(dat_30, subgroup %in% cond4)
dat_30_goi <- dat_30_goi[c(-4, -11, -16, -21), ]
dat_30_goi <- dat_30_goi[-16, ] # if include spc group
dat_30_goi <- dat_30_goi[-11, ]
log.dat.30 <- log(dat_30_goi[,3:14])






groups <- dat_30_goi[,1]
subgroups <- dat_30_goi[,2]


pca <- prcomp(log.dat.30,
              center = T,
              scale. = T)
plot(pca, type = 'l')
summary(pca)

#11 sample in first-order is lame, also 4 if mixing with spc
#16 sample in spc is lame

g <- ggbiplot(pca, obs.scale = 1, var.scale = 1, groups=subgroups, ellipse = T, circle = F)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')

print(g)

dataPCA <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], condition = subgroups)
z <- qplot(PC1, PC2, data = dataPCA, color = condition)
z + geom_label(aes(label = rownames(dat_30_goi)))


#90 minutes data
dat_90_goi <- filter(dat_90, subgroup %in% cond1)
dat_90_goi <- dat_90_goi[-15,]
log.dat.90 <- log(dat_90_goi[,3:14])

#lame 15 sample

subgroups <- dat_90_goi[,2]

pca <- prcomp(log.dat.90,
              center = T,
              scale. = T)
plot(pca, type = 'l')
summary(pca)

subgroups <- factor(c(rep("2nd order", 7), rep("Control", 7)), levels = c("2nd order", "Control"))

tiff(file = "/home/arthur/Documents/R_learning/biplot.tiff", width = 3500, height = 4400, units = "px", res = 600)
pc <- princomp(log.dat.90, cor = T, score = T)
g <- ggbiplot(pca, obs.scale = 1, var.scale = 1, groups=subgroups, ellipse = T, circle = F, varname.size = 4)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', legend.position = 'top', legend.text = element_text(size = 16), axis.title = element_text(size = 18), axis.text = element_text(size = 18))
g <- g + geom_vline(xintercept = 0, linetype = 'dashed', color = 'red')
print(g)
dev.off() 

loadings(pc)

?ggbiplot

#270 minutes data

dat_270_goi <- filter(dat_270, subgroup %in% cond1)
log.dat.270 <- log(dat_270_goi[,3:14])


subgroups <- dat_270_goi[,2]


pca <- prcomp(log.dat.270,
              center = T,
              scale. = T)


plot(pca, type = 'l')
summary(pca)

pc <- princomp(log.dat.270, cor = T, score = T)
g <- ggbiplot(pca, obs.scale = 1, var.scale = 1, groups=subgroups, ellipse = F, circle = F)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')

print(g)

pc$loadings

# Correlations 

res <- cor(dat_270[,c(3:14)])
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = res, col = col, symm = TRUE, na.rm = T)


# The whole dataset
dat_quant <- na.omit(dat)

sec_ord <- dat_quant[dat_quant$subgroup == "second-order",]
res_sec_ord <- cor(sec_ord[,-c(1,2)])
chart.Correlation(res_sec_ord, histogram=TRUE, pch=19)

fir_ord <- dat_quant[dat_quant$subgroup == "first-order",]
res_fir_ord <- cor(fir_ord[,-c(1,2)])
chart.Correlation(res_fir_ord, histogram=TRUE, pch=19)

spc <- dat_quant[dat_quant$subgroup == "spc",]
res_spc <- cor(spc[,-c(1,2)])
chart.Correlation(spc, histogram=TRUE, pch=19)

heatmap(x = res_spc, col = col, symm = TRUE, na.rm = T)


col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = res_sec_ord, col = col, symm = TRUE, na.rm = T)
heatmap(x = res_fir_ord, col = col, symm = TRUE, na.rm = T)



turned_fir_ord <- t(dat_30[,-c(1,2)])
nam <- dat_30$subgroup
names(turned_fir_ord) <- nam
turned_fir_ord <- as.data.frame(turned_fir_ord)
distance_fo <- dist(turned_fir_ord)
hc.c.fo <- hclust(distance_fo)
plot(hc.c.fo, labels = rownames(turned_fir_ord))


turned_sec_ord <- t(dat_90[,-c(1,2)])
nam <- dat_90$subgroup
names(turned_sec_ord) <- nam
turned_sec_ord <- as.data.frame(turned_sec_ord)
distance_so <- dist(turned_sec_ord)
hc.c.so <- hclust(distance_so)
plot(hc.c.so, labels = rownames(turned_sec_ord))





spc <- dat_quant[c(16:23),]
spc <- dat_quant[c(1:8),]
res <- cor(spc)
chart.Correlation(spc, histogram=TRUE, pch=19)


# distance <- dist(dat)
distance







my_data <- mtcars[, c(1,3,4,5,6,7)]

chart.Correlation(res, histogram=TRUE, pch=19)


res <- cor(dat_30[,c(3:14)])
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = res, col = col, symm = TRUE, na.rm = T)

res <- cor(dat_90[,c(3:14)])
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = res, col = col, symm = TRUE, na.rm = T)

round(res, 2)
res








dat_90_goi <- dat_90[dat_30$subgroup == "second-order" | dat_30$subgroup == "first-order", ]
dat_90_goi <- dat_90_goi[-15,]
log.dat.90 <- log(dat_90_goi[,3:14])


subgroups <- dat_90_goi[,2]

pca <- prcomp(log.dat.90,
              center = T,
              scale. = T)
plot(pca, type = 'l')
summary(pca)

pc <- princomp(log.dat.270, cor = T, score = T)

library(ggbiplot)

g <- ggbiplot(pca, obs.scale = 1, var.scale = 1, groups=subgroups, ellipse = T, circle = F)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')

print(g)




dat_goi <- dat[dat$subgroup == "second-order" | dat$subgroup == "first-order", ]
log.dat <- log(dat_goi[,3:13])


subgroups <- dat_goi[,2]

pca <- prcomp(log.dat,
              center = T,
              scale. = T)
plot(pca, type = 'l')
summary(pca)


g <- ggbiplot(pca, obs.scale = 1, var.scale = 1, groups=subgroups, ellipse = T, circle = F)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')

print(g)

dataPCA <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], condition = subgroups)
z <- qplot(PC1, PC2, data = dataPCA, color = condition)
plot(dataPCA$PC1, dataPCA$PC2)
text( pca$x[,1], pca$x[,2], rownames( dat ), pos= 3 )

z + geom_label(aes(label = rownames(dat_goi)))
z


# http://forrest.psych.unc.edu/research/vista-frames/help/lecturenotes/lecture13/biplot.html   Explanation of PCA fpr dummies



var <- get_pca_var(pca)
var$contrib
fviz_pca_var(pca, col.var = "black")
fviz_contrib(pca, choice = "var", axes = 1:2, top = 12)


tiff(file = "/home/arthur/Documents/R_learning/contr.tiff", width = 4000, height = 4400, units = "px", res = 600)
contr_plot <- fviz_pca_var(pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

contr_plot + theme( plot.title = element_text(size = 18), legend.title = element_text(size = 18), legend.text = element_text(size = 14), axis.title = element_text(size = 18), axis.text = element_text(size = 18))
dev.off()

?princomp
loadings(pca)




summary(pc.cr <- princomp(USArrests, cor = TRUE))
loadings(pc.cr)  # note that blank entries are small but not zero
## The signs of the columns are arbitrary
plot(pc.cr) # shows a screeplot.
biplot(pc.cr)

summary(pc.cr <- princomp(dat_30_goi[,-1], cor=TRUE))



library(cluster)
library(HSAUR)
data(pottery)

km    <- kmeans(turned_fir_ord,3)
dissE <- daisy(turned_fir_ord) 
dE2   <- dissE^2
sk2   <- silhouette(km$cl, dE2)
plot(sk2)

clusplot(turned_fir_ord, km$cluster, color=TRUE, shade=T,   lines=0)





# PRC summary
setwd("/home/arthur/Documents/R_learning")
dat_30 <- read.xlsx("for PCA.xlsx", sheetIndex = 1, header = T)
dat_90 <- read.xlsx("for PCA.xlsx", sheetIndex = 2, header = T)
dat_270 <- read.xlsx("for PCA.xlsx", sheetIndex = 3, header = T)

cond2 <- c("second-order", 'first-order', "spc")
cond1 <- c("Second-order", 'First-order')
cond3 <- c('first-order', "spc")
cond4 <- c("Second-order", "spc")


pca_table <- function(dataset, condit){
  groups <- dataset[,1]
  dat_goi <- filter(dataset, groups %in% condit)
  log.dat <- log(dat_goi[,2:14])
  groups <- dat_goi[,1]
  pca <- prcomp(log.dat[,-1],
                center = T,
                scale. = T)
  pca_load <- princomp(log.dat[,-1], cor = TRUE)
  newList <- list("pca" = pca, "groups" = groups, "load" = pca_load)
  return(newList)
  
}


pca_plot <- function(pca){
  g <- ggbiplot(pca$pca, obs.scale = 1, var.scale = 1, groups=pca$groups, ellipse = T, circle = F)
  g <- g + scale_color_discrete(name = '')
  g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
  print(g)
}
  
pca_table(dat_30, cond4)  


pca_contribution <- function(pca){
  var <- get_pca_var(pca)
  contr <- var$contrib
  print(contr)
  contr_plot <- fviz_pca_var(pca, col.var = "contrib",
                            gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
  
  contr_plot + theme( plot.title = element_text(size = 18), legend.title = element_text(size = 18), 
                      legend.text = element_text(size = 14), axis.title = element_text(size = 18), 
                      axis.text = element_text(size = 18))
}

pca <- pca_table(dat_270, cond4)
pca_plot(pca)
pca_contribution(pca$pca)
loadings(pca$load)
