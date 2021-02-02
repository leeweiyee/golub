# Lab 5: Cluster Analysis Using R and Bioconductor
# Weiyee Lee
# 19 Jan 2021

library(Biobase)
library(annotate)
library(golubEsets)
library(genefilter)
library(ellipse)
library(lattice)
library(cluster)
library(MVA)

# get train data
data(Golub_Train)
# extract expression matrix
X <- exprs(Golub_Train)
X[X < 100] <- 100
X[X > 16000] <- 16000

# filter
mmfilt <- function(r = 5, d = 500, na.rm = TRUE) {
  function(x) { 
    minval <- min(x, na.rm = na.rm) 
    maxval <- max(x, na.rm = na.rm)
    (maxval/minval > r) && (maxval - minval > d)
  }
}
mmfun <- mmfilt()
ffun <- filterfun(mmfun) #?
sub <- genefilter(X, ffun)
sum(sub) # 3051

# subset
X <- X[sub, ]
X <- log2(X)
golubTrainSub <- Golub_Train[sub, ]
exprs(golubTrainSub) <- X #?
Y <- golubTrainSub$ALL.AML
Y <- paste(Golub_Train$ALL.AML, Golub_Train$T.B.cell)
Y <- sub("NA", "", Y)

# correlation distance
r <- cor(X)
dimnames(r) <- list(as.vector(Y), as.vector(Y))
d <- 1 - r

# plot correlation matrix 1
plotcorr(r, main="Leukemia data: Correlation matrix for 38 mRNA samples\n All 3,051 genes")

# plot correlation matrix 2
plotcorr(r,numbers=TRUE, main="Leukemia data: Correlation matrix for 38 mRNA samples\n All 3,051 genes")

# plot correlation matrix 3
levelplot(r,col.region=heat.colors(50), main="Leukemia data: Correlation matrix for 38 mRNA samples\n All 3,051 genes")

# selecting genes where p-value < 0.01
gtt <- ttest(golubTrainSub$ALL, p = 0.01)
gf1 <- filterfun(gtt)
whT <- genefilter(golubTrainSub, gf1)
sum(whT) # 609
# subset genes
gTrT <- golubTrainSub[whT, ]

# plot correlation of subset genes 1
rS <- cor(exprs(gTrT))
dimnames(rS)<-list(gTrT$ALL, gTrT$ALL)
dS <- 1 - rS
levelplot(rS,col.region=heat.colors(50),main="Leukemia data: Correlation matrix for 38 mRNA samples\n 609 genes")

# selecting genes where p-value < 0.05/3051
gtt <- ttest(golubTrainSub$ALL, p = 0.05/3051)
gf1 <- filterfun(gtt)
whT <- genefilter(golubTrainSub, gf1)
sum(whT) # 74
# subset genes
gTrT <- golubTrainSub[whT, ]

# plot correlation of subset genes 2
rS <- cor(exprs(gTrT))
dimnames(rS)<-list(gTrT$ALL, gTrT$ALL)
dS <- 1 - rS
levelplot(rS,col.region=heat.colors(50),main="Leukemia data: Correlation matrix for 38 mRNA samples\n 74 genes")

# multidimensional scaling (MDS) 1
mds <- cmdscale(d, k = 2, eig = TRUE)
plot(mds$points, type = "n", xlab = "", ylab = "", main = "MDS for ALL AML data, correlation matrix, G=3,051 genes, k=2")
text(mds$points[, 1], mds$points[, 2], Y, col = as.integer(factor(Y)) + 1, cex = 0.8)

# multidimensional scaling (MDS) 2
mds3 <- cmdscale(d, k = 3, eig=TRUE)
pairs(mds3$points, main="MDS for ALL AML data, correlation matrix, G=3,051 genes, k=3", pch=c("B","T","M")[as.integer(factor(Y))], col = as.integer(factor(Y))+1)

# does not work
if( require(rgl) && interactive() ) { 
  rgl.spheres(x=mds3$points[,1], mds3$points[,2], mds3$points[,3], col= as.integer(factor(Y))+1, radius=.01) }

# scree plot
mdsScree <- cmdscale(d, k = 8, eig = TRUE)
plot(mdsScree$eig, pch = 18, col = "blue")

# expression density diagnostics (EDD)
data(Golub_Merge)
gMs <- Golub_Merge[,(Golub_Merge$T.B != "T-cell" | is.na(Golub_Merge$T.B))]
fF1 <- gapFilter(100, 200, .1)
ff <- filterfun(fF1)
whgMs <- genefilter(gMs, ff)
sum(whgMs) ##3082
gMs <- gMs[whgMs,]

library(edd)
gMsALL <- gMs[,gMs$ALL=="ALL"]
gMsAML <- gMs[,gMs$ALL=="AML"]
#set the seed for the rng for reproducibility
set.seed(1234)
edd.ALL <- edd(gMsALL, method="knn", k=4, l=2) # doubt
sum(is.na(edd.ALL)) # 29
edd.AML <- edd(gMsAML, method="knn", k=4, l=2)
sum(is.na(edd.AML)) # 71
table(edd.ALL, edd.AML)

# geneComp
diff1 <- edd.AML == "N(0,1)" & edd.ALL == ".75N(0,1)+.25N(4,1)"
table(diff1)
diff1[is.na(diff1)] <- FALSE 
AMLexprs <- exprs(gMsAML)[diff1,]
ALLexprs <- exprs(gMsALL)[diff1,]

# transform
tAML <- fq.matrows(t(apply(AMLexprs, 1, centerScale)))
tALL <- fq.matrows(t(apply(ALLexprs, 1, centerScale)))
par(mfrow=c(2,2))
hist(AMLexprs[1,])
hist(ALLexprs[1,])
hist(AMLexprs[2,])
hist(ALLexprs[2,])
par(mfrow=c(1,1)) # reset plotting parameter

# agglomerative hierarchical clustering: average linkage
hc1 <- hclust(as.dist(d), method="average")
coph1 <- cor(cophenetic(hc1),as.dist(d))
plot(hc1, main=paste("Dendrogram for ALL AML data: Coph = ", round(coph1,2)), sub="Average linkage, correlation matrix, G=3,051 genes")
cthc1 <- cutree(hc1, 3)
table(Y, cthc1)

# agglomerative hierarchical clustering: single linkage
hc2 <- hclust(as.dist(d), method="single")
coph2 <- cor(cophenetic(hc2),as.dist(d))
plot(hc2, main=paste("Dendrogram for ALL AML data: Coph = ", round(coph2,2)), sub="Single linkage, correlation matrix, G= 3,051 genes")
cthc2 <- cutree(hc2, 3)
table(Y, cthc2)

# agglomerative hierarchical clustering: complete linkage
hc3 <- hclust(as.dist(d), method="complete")
coph3 <- cor(cophenetic(hc3),as.dist(d))
plot(hc3, main=paste("Dendrogram for ALL AML data: Coph = ", round(coph3,2)), sub="Complete linkage, correlation matrix, G= 3,051 genes")
cthc3 <- cutree(hc3, 3)
table(Y, cthc3)

# divisive hierarchical clustering
di1 <- diana(as.dist(d))
cophdi <- cor(cophenetic(di1), as.dist(d))
plot(di1, which.plots=2, #$ main=paste("Dendrogram for ALL AML data: Coph = ", round(cophdi,2)), sub="Divisive algorithm, correlation matrix, G= 3,051 genes")
ct.di <- cutree(di1, 3)
table(Y, ct.di)

# partitioning around medoids (PAM) with three clusters
set.seed(12345)
pm3 <- pam(as.dist(d), k=3,diss=TRUE)
table(Y, pm3$clustering)
# cluster plot
clusplot(d, pm3$clustering, diss=TRUE, labels=3, col.p=1, col.txt=as.integer(factor(Y))+1, main="Bivariate cluster plot for ALL AML data\n Correlation matrix, K=3, G=3,051 genes")
plot(pm3,which.plots=2, main="Silhouette plot for ALL-AML Data")

# PAM with four clusters
set.seed(12345)
pm4 <- pam(as.dist(d), k=4, diss=TRUE)
clusplot(d, pm4$clustering, diss=TRUE, labels=3, col.p=1, col.txt=as.integer(factor(Y))+1, main="Bivariate cluster plot for ALL AML data\n Correlation matrix, K=4, G=3,051 genes")

