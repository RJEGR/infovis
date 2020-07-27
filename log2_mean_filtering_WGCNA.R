library(DESeq2)

path <- c('~/transcriptomics/Diana_Lara/Diana_results/')

count_arg <- "iso.counts.matrix$"

MIN_CPM <- 1

count_file <- list.files(path = path, pattern = count_arg, 
                         full.names = T)

count <- read.table(count_file, header=T, com='', 
                    row.names=1, check.names=F, sep='\t', 
                    stringsAsFactors = FALSE)

# mean_threshold <- round(mean(rowSums(edgeR::cpm(count))))
mean_threshold <- round(mean(rowSums(count)))

MIN_CPM <- mean_threshold

qntl <- c(0., 0.25, 0.5, 0.75, 0.99, 1.0)

prep_boxdata <- function(x) as.numeric(quantile(x, qntl, na.rm=T))
boxdata <- log2(apply(count, 2, prep_boxdata) + 1)
boxplot(boxdata)
abline(h = log2(MIN_CPM + 1), col = 'red')

count <- count[rowSums(count > MIN_CPM) >= MIN_REPS, ]

conditions <- sapply(strsplit(names(count), "_"), `[`, 2)

conditions <- data.frame(conditions=factor(conditions))

rownames(conditions) <- colnames(count)

count <- round(count)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = count,
  colData = conditions,
  design = ~ conditions )

dds <- estimateSizeFactors(ddsFullCountTable)

datExpr <- counts(dds, normalized = TRUE)
datExpr <- t(log2(datExpr + 1))

saveRDS(datExpr, file = paste0(path, '/', 'datExpr.rds'))
saveRDS(count, file = paste0(path, '/', 'count_prepared.rds'))

