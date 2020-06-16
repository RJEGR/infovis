library(DESeq2)
library(edgeR)

genes <- 1*10^4
x <- DESeq2::makeExampleDESeqDataSet(n= genes, m = 10)

colData <- colData(x)

gene.length <- sample(rep(100:genes), genes, replace = T)


# L <- mean(x$sizeFactor) * 1e-6
# M <- median(x$sizeFactor) * 1e-6
# c(L, M)
# 
# lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(x)

samplenames <- rownames(colData)

col <- brewer.pal(nsamples, "Paired")

par(mfrow=c(1,3))

lcount <- log(DESeq2::counts(x))
plot(density(lcount[,1]), col=col[1], lwd=2, 
     ylim=c(0,1), 
     las=2, main="", xlab="")
title(main="A. Raw data", xlab="")
# abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcount[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

x2 <- estimateSizeFactors(x)
lcpm <- cpm(counts(x2), log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, 
     ylim=c(0,1), 
     las=2, main="", xlab="")
title(main="B. cpm", xlab="Log-scale")
# abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

# 3
lrpkm <- rpkm(counts(x2), gene.length = gene.length,
                 normalized.lib.sizes = TRUE,
                 prior.count = 2, log = T)

plot(density(lrpkm[,1]), col=col[1], lwd=2, 
     ylim=c(0,1), 
     las=2, main="", xlab="")
title(main="C. RPKM", xlab="")
#abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lrpkm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

# normalization ----

par(mfrow=c(1,2))
lcpm <- cpm(counts(x), log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Unnormalised data",ylab="Log-scale")
# 
# x2 <- estimateSizeFactors(x)
# lcpm <- cpm(counts(x3), log=TRUE)
# boxplot(lcpm, las=2, col=col, main="")
# title(main="B. Normalised data",ylab="Log-cpm")
x2 <- estimateSizeFactors(x)
lrpkm <- rpkm(counts(x2), gene.length = gene.length,
              normalized.lib.sizes = TRUE,
              prior.count = 2, log = T)
boxplot(lrpkm, las=2, col=col, main="")
title(main="C. Normalised rpkm",ylab="Log-scale")

# pca

counts <- DESeq2::counts(x)

mtd <- cbind(samplenames, colData)

PCA <- prcomp(t(lrpkm), scale. = FALSE)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dtvis <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                    mtd)

pca_plot2 <- ggplot(dtvis, aes(PC1, PC2)) +
  #geom_point(aes(shape = condition), size = 5, alpha = 0.9) +
  geom_label(aes(label = condition, fill = samplenames), alpha = 0.9) +
  labs(caption = '') +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  #coord_fixed(ratio = sd_ratio) +
  scale_color_manual(name = 'Samples', values = col)  +
  scale_fill_manual(name = '', values = col)  +
  scale_shape_discrete(name = '') +
  guides(color = FALSE) +
  #scale_color_manual(name=NULL, values=color_sample) +
  theme_bw()

library(patchwork)

pca_plot + pca_plot2
#

# sample corr
library(corrplot)


par(mfrow=c(1,2))
# corrplot(cor(lcpm), title = '')
corrplot(cor(lrpkm), title = 'Normalised rpkm', type = 'upper')
corrplot(cor(counts(x)), title = 'Unnormalised-data', type = 'upper')


