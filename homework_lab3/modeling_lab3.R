library(devtools)
library(Biobase)
library(limma)
library(edgeR)
library(genefilter)
library(qvalue)
library(tidyverse)
library(data.table)

load(file="bottomly.Rdata")
ls()

edata <- as.matrix(exprs(bottomly.eset))
dim(edata)
edata[1:5,1:5]

edata <- log2(as.matrix(edata) + 1)
edata <- edata[rowMeans(edata) > 10, ]

library(RColorBrewer)
library(gplots)

my_palette <- colorRampPalette(c("blue", "white", "orange"))(n = 299)

png("bottomly_heatmap_raw.png",height=700,width=700)
heatmap.2(edata,
          main = "Bottomly et al. Raw", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          dendrogram="none",     # only draw a row dendrogram
          scale = "row",
          Colv=FALSE)
dev.off()


png("bottomly_heatmap_clustered.png",height=700,width=700)
heatmap.2(edata,
          main = "Bottomly et al. Clustered", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          dendrogram="none",     # only draw a row dendrogram
          scale = "row")
dev.off()

# problem1
pdf("gozdziewska_problem1.pdf",height=1000,width=1000)
heatmap.2(edata,
          main = "Bottomly et al. Clustered with column dendrogram", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none", # turns off trace lines inside the heat map
          key=T,
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          dendrogram="col",     # only draw a row dendrogram
          scale = "column") 
dev.off()


# SVD
edata <- t(scale(t(edata), scale=FALSE, center=TRUE))
svd.out <- svd(edata)
names(svd.out)

print(paste("Dimension of left singular vectors:", dim(svd.out$u)))
print(paste("Length of singular values:",length(svd.out$d)))
print(paste("Dimension of right singular vectors:",dim(svd.out$v)))

par(mfrow=c(1,2))
plot(svd.out$d, pch=20, ylab="Singular values")
plot(svd.out$d^2/sum(svd.out$d^2)*100, pch=20, ylab="% variance explained")

plot(1:ncol(edata), svd.out$v[,1],pch=20)
plot(1:ncol(edata), svd.out$v[,2],pch=20)
plot(1:ncol(edata), svd.out$v[,3],pch=20)

PC = data.table(svd.out$v,pData(bottomly.eset))
ggplot(PC) + geom_point(aes(x=V1, y=V2, col=as.factor(strain)))
ggplot(PC) + geom_point(aes(x=V1, y=V2, col=as.factor(lane.number)))

ggplot(PC) + geom_point(aes(x=V1, y=V2, col=as.factor(experiment.number)))

# problem2
ggplot(PC) + geom_point(aes(x=V3, y=V4, col=as.factor(strain)))

ggplot(PC) + geom_boxplot(aes(x=as.factor(strain), y=V1))
ggplot(PC) + geom_violin(aes(x=as.factor(strain), y=V1),draw_quantiles = c(0.25, 0.5, 0.75))
ggplot(PC) + geom_violin(aes(x=as.factor(strain), y=V1),draw_quantiles = c(0.25, 0.5, 0.75)) + geom_jitter(aes(x=as.factor(strain), y=V1))

# problem3
PC = data.table(svd.out$u,pData(bottomly.eset))
ggplot(PC) + geom_point(aes(x=V1, y=V2, col=as.factor(strain)))

# problem4
PC = data.table(svd.out$u,pData(bottomly.eset))[1:5]
ggplot(PC) + geom_violin(aes(x=as.factor(strain), y=V1),draw_quantiles = c(0.25, 0.5, 0.75)) + geom_jitter(aes(x=as.factor(strain), y=V1))

pc1 = prcomp(edata)
plot(pc1$rotation[,1],svd.out$v[,1])

edata.col <- scale(edata, scale=FALSE, center=TRUE)
svd.col <- svd(edata.col)
plot(pc1$rotation[,1],svd.col$v[,1],col=2)
abline(0,1)

all(pc1$rotation[,1] == svd.col$v[,1])

library(irlba)
tsvd.out <- irlba(edata, nv = 4)
dim(tsvd.out$u)
length(tsvd.out$d)
dim(tsvd.out$v)
plot(tsvd.out$v[,1],-1*svd.out$v[,1]); abline(0,1,col="red")
plot(tsvd.out$v[,2],svd.out$v[,2]); abline(0,1,col="red")


library(irlba)
library(Rtsne)

set.seed(1)
tsne_out <- Rtsne(edata,pca=TRUE,perplexity=30)
tsne_out = data.table(tsne_out$Y)
ggplot(tsne_out) + geom_point(aes(x=V1, y=V2))

# problem5
cluster <- kmeans(edata, centers = 5)
set.seed(1)
tsne_out <- Rtsne(edata,pca=TRUE,perplexity=30)
tsne_out = data.table(tsne_out$Y)
ggplot(data=tsne_out) + geom_point(aes(x=V1, y=V2, color=factor(cluster$cluster)))
                              