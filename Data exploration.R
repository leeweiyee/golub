# Golub dataset exploration
# Weiyee Lee
# 25 Jan 2021

library(Biobase)
library(annotate)
library(golubEsets)
library(genefilter)
library(ellipse)
library(lattice)
library(cluster)
library(MVA)

class(Golub_Train) # ExpressionSet class
data <- exprs(Golub_Train)
nrow(data) # 7129 genes
ncol(data) # 38 samples
min(data) # -28400
max(data) # 61228
data[data < 100] <- 100
data[data > 16000] <- 16000

# filter genes
mmfilt <- function(r = 5, d = 500, na.rm = TRUE) {
  function(x) { 
    minval <- min(x, na.rm = na.rm) 
    maxval <- max(x, na.rm = na.rm)
    (maxval/minval > r) && (maxval - minval > d)
  }
}
mmfun <- mmfilt()
ffun <- filterfun(mmfun) #?
sub <- genefilter(data, ffun)
sum(sub) # 3051

# get subset
data <- data[sub, ]
nrow(data) # 3051 filtered genes
data <- log2(data) # normalise data
min(data) # 6.643856
max(data) # 13.96578
golubTrainSub <- Golub_Train[sub, ] # sub-expressionset
exprs(golubTrainSub) <- data # seems redundant?

labels <- golubTrainSub$ALL.AML
labels <- paste(Golub_Train$ALL.AML, Golub_Train$T.B.cell)
labels <- sub("NA", "", labels) # edit sample labels

library(pheatmap)
library(RColorBrewer)
library(viridisLite)

# 74 differentially expressed (DE) genes
gtt <- ttest(golubTrainSub$ALL, p = 0.05/3051)
gf1 <- filterfun(gtt)
whT <- genefilter(golubTrainSub, gf1)
sum(whT) # 74 genes
gTrT <- golubTrainSub[whT, ]
a1 <- exprs(gTrT)

# 3051-74=2977 non-DE genes
whT <- genefilter(golubTrainSub, !gf1)
sum(!whT) # 2977 genes
gTrT <- golubTrainSub[!whT, ]
b1 <- exprs(gTrT)

nrow(a1) # 74 DE genes
nrow(b1) # 2977 non-DE genes

# heatmap of 74 DE genes
sample_columns <- data.frame('Tumour Type' = labels)
rownames(sample_columns) <- colnames(a)
mat_colors <- list(group = brewer.pal(3, "Set1"))
names(mat_colors$group) <- unique(labels)
pheatmap(mat = a1,
         border_color = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_col = sample_columns,
         annotation_colors = mat_colors,
         drop_levels = TRUE,
         fontsize = 9,
         main = "Leukemia data: Heatmap for 38 mRNA samples\n 74 genes")

# 609 differentially expressed (DE) genes
gtt <- ttest(golubTrainSub$ALL, p = 0.01)
gf1 <- filterfun(gtt)
whT <- genefilter(golubTrainSub, gf1)
sum(whT) # 609 genes
gTrT <- golubTrainSub[whT, ]
a2 <- exprs(gTrT)

# 3051-609=2442 non-DE genes
gtt <- ttest(golubTrainSub$ALL, p = 0.01)
gf1 <- filterfun(gtt)
whT <- genefilter(golubTrainSub, !gf1)
sum(!whT) # 2442 genes
gTrT <- golubTrainSub[!whT, ]
b2 <- exprs(gTrT)

nrow(a2) # 74 DE genes
nrow(b2) # 2977 non-DE genes

# heatmap of 609 differentially expressed (DE) genes
sample_columns <- data.frame('Tumour Type' = labels)
rownames(sample_columns) <- colnames(a)
mat_colors <- list(group = brewer.pal(3, "Set1"))
names(mat_colors$group) <- unique(labels)
pheatmap(mat = a2,
         border_color = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_col = sample_columns,
         annotation_colors = mat_colors,
         drop_levels = TRUE,
         fontsize = 9,
         main = "Leukemia data: Heatmap for 38 mRNA samples\n 609 genes")

# box plot of random DE and non-DE genes
rand_genes_de <- sample.int(nrow(a1), 2)
rand_gene_1 <- a1[rand_genes_de[1], ]
rand_gene_2 <- a1[rand_genes_de[2], ]

rand_genes_nde <- sample.int(nrow(b1), 2)
rand_gene_3 <- b1[rand_genes_nde[1], ]
rand_gene_4 <- b1[rand_genes_nde[2], ]

rand_gene_1df <- data.frame(expression = rand_gene_1,
                            tumour_type = labels)
rand_gene_2df <- data.frame(expression = rand_gene_2,
                            tumour_type = labels)
rand_gene_3df <- data.frame(expression = rand_gene_3,
                            tumour_type = labels)
rand_gene_4df <- data.frame(expression = rand_gene_4,
                            tumour_type = labels)
library(ggplot2)

box1 <- ggplot(rand_gene_1df, 
               aes(x = tumour_type, y = expression, fill = tumour_type)) +
  theme(axis.title.x=element_blank()) +
  labs(fill = NULL) +
  geom_boxplot()
box2 <- ggplot(rand_gene_2df, 
               aes(x = tumour_type, y = expression, fill = tumour_type)) + 
  theme(axis.title.x=element_blank()) +
  labs(fill = NULL) +
  geom_boxplot()
box3 <- ggplot(rand_gene_3df, 
               aes(x = tumour_type, y = expression)) + 
  theme(axis.title.x=element_blank()) +
  labs(fill = NULL) +
  geom_boxplot()
box4 <- ggplot(rand_gene_4df, 
               aes(x = tumour_type, y = expression)) + 
  theme(axis.title.x=element_blank()) +
  labs(fill = NULL) +
  geom_boxplot()

library(plotly)

subplot(box1, box2, box3, box4, 
        nrows = 2) %>%
  layout(title = 'Leukemia data: expression of 2 random DE and 2 random non-DE genes')

# scatter plot of random 74 DE genes
rand_genes_de <- sample.int(nrow(a1), 2)
rand_gene_1 <- a1[rand_genes_de[1], ]
rand_gene_2 <- a1[rand_genes_de[2], ]

plot(rand_gene_1, rand_gene_2, main = 'Leukemia data: expression of 2 random DE genes')
abline(reg = lm(rand_gene_2 ~ rand_gene_1))

# PCR
library(ggfortify)

pca_res <- prcomp(t(data[whT,]), scale. = TRUE)
autoplot(pca_res)
df <- as.data.frame(cbind(t(data[whT,]), labels))
autoplot(pca_res, data = df, colour = 'labels')
