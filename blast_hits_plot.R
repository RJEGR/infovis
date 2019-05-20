# 1 meta-assembly of oyster samples tissues
# 2 transform to super-transcript each isoforms-gene group
# 3  predict ORF with transdecoder

# 4 BLASTP

# query=trinity_genes.fasta.transdecoder.pep
# reference=peerj-04-1763-s005.pep
#blastp -query $query \
#-db $DBS/$reference -num_threads 24 \
#-max_target_seqs 1 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen" -evalue 1e-5 \
#> $filename


rm(list=ls()); 
options(stringsAsFactors = FALSE)

# Load package ----
.cran_packages <- c('dplyr', 'superheat', 'ggplot2', 'RColorBrewer')
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst], dep=TRUE, repos='http://cran.us.r-project.org')
}
sapply(c(.cran_packages), require, character.only = TRUE)

# Reading
# BLASTP machinery in transcriptome ----

blast_load <- function(x) { 
  blast_obj <- read.table(x, stringsAsFactors = FALSE)
  blast_names <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs", "qcovhsp", "qlen", "slen")
  names(blast_obj) <- blast_names
  #blast_obj <- blast_obj[, c("pident", "mismatch","gapopen", "qcovs", "length")]
  return(blast_obj)
}


blast_path = '/Users/cigom/transcriptomics/miRNAs/blastp'
blast_files <- list.files(blast_path, full.names = TRUE, pattern = 'trinity_genes.fasta.transdecoder_vs_peerj-04-1763-s005_blastp.outfmt6')

blast_plot <- do.call(rbind,lapply(blast_files, blast_load))
blast_plot$qseqid <- sapply(strsplit(blast_plot$qseqid, "[::]"), `[`, 3)

blast_plot$sseqid <- make.unique(as.vector(blast_plot$sseqid), sep = "_")

#blastn_plot$name <- factor(blastn_plot$name, levels = c('1e-10', '1e-20', '1e-40', '1e-80', '1e-120', '1e-200', '1e-300'))

# Subset blastp ----

n_sseqid <- 10

ssqid_text <- names(table(blast_plot$sseqid)[table(blast_plot$sseqid) > n_sseqid])

# Viz blastp ----
path <- "/Users/cigom/transcriptomics/miRNAs"
png(paste0(path,"/", "aligment_profile.png"), units="px", width=2400, height=2600, res=400)

# Change point shapes, colors and sizes, shape=cyl,
ggplot(blast_plot, aes(x=pident, y=length)) +
  geom_point(aes(color = qcovs, size = mismatch), alpha = 0.7) +
  # scale_shape_manual(values = 1:length(levels(blastn_plot$name))) + # if shape = name
  scale_color_distiller(palette = "Purples", direction = 1) +
  theme_classic() +
  labs(x = 'Identidad del alineamiento (%)',
       y = 'Longitud del alineamiento',
       subtitle = 'Blast hits (blastp): Oyster pseudo-genes vs miRNA machinery DataBase',
       caption = paste0('miRNA machinery with Identity >= 90 & Coverage >= 90 is labeled')) +
       #caption = paste0('miRNA machinery with Frequency > ', n_sseqid, " is labeled")) +
  geom_text(data = subset(blast_plot, pident >= 90 & qcovs >= 90 ), aes(label=sseqid), hjust = 1, vjust = 1, check_overlap =TRUE, size = 2.3)
  #geom_text(data = subset(blast_plot, sseqid==ssqid_text), aes(label=sseqid), hjust = 0, vjust = 1, check_overlap =TRUE, size = 2.3)
dev.off()
  #facet_wrap( ~ name, scales = 'free')

# ggplot() +
#   stat_summary_2d(data = blast_plot,
#                   mapping = aes(x = qcovs,
#                                 y = pident,
#                                 z = mismatch),
#                   fun = mean)
# Output subset list ----

write(subset(blast_plot, pident >= 90 & qcovs >= 90 )$qseqid,
      file = paste0(blast_path, '/', 'blastp_pident_90_qcovs_90.list'))

# grep -A 3 -Ff blastp_pident_90_qcovs_90.list trinity_genes.fasta.transdecoder.cds | sed '/^--$/d' | awk '/^>Gene/{gsub(/[::]/, " "); print ">"$2; next}{print}' > trinity_genes.fasta.transdecoder.pident_90_qcovs_90.cds
# grep "^>" trinity_genes.fasta.transdecoder.pident_90_qcovs_90.cds | sed 's/>//g' > blastp_pident_90_qcovs_90.genes.ids

# Parse to RMPK ----
#

count_path <- "~/transcriptomics/miRNAs/subset_abundance/"
count_file <- list.files(count_path, full.names = TRUE, pattern = "RSEM.isoforms.gene.counts.matrix")
count_matrix <- readr::read_delim(count_file, delim = "\t")


samples_file <- list.files(count_path, full.names = TRUE, pattern = "MetadataCv_exposure.csv")
samples_data <- read.table(samples_file, sep = ',', header = TRUE)
samples_data = samples_data[samples_data[,2] != '',]
names(samples_data)[1] <- c("sample_name")
# es una excelent idea transformar toda la matriz de datos y excluir expresiones menores a n % 
## https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#data-transformations-and-visualization
# library("DESeq2")
# cts <- as.matrix(count_matrix[-1])
# 
# rownames(cts) <- count_matrix$id
# 
# samples_file <- list.files(count_path, full.names = TRUE, pattern = "samples.file")
# coldata <- read.table(samples_file)
# 
# rownames(coldata) <- coldata$V2
# colnames(coldata)[1] <- "condition"
# all(rownames(coldata) == colnames(cts))
# cts <- cts[, rownames(coldata)]
# all(rownames(coldata) == colnames(cts))
# dds <- DESeqDataSetFromMatrix(countData = cts,
#                               colData = coldata,
#                               design = ~ condition)
# dds


names(count_matrix)[1] <- 'id'

# blast_subset <- blast_plot[blast_plot$pident >= 90 & blast_plot$qcovs >= 90, ]

# x <- blast_plot[blast_plot$pident >= 90 & blast_plot$qcovs >= 90, c('qseqid', 'sseqid')]

x <- blast_plot

blast_subset_ids <- data.frame(id = x[,1], name = x[,2],
                               pident = x$pident, 
                               qcovs = x$qcovs,
                               stringsAsFactors = FALSE)

count_matrix %>% 
  group_by(id) %>%
  inner_join(blast_subset_ids, by = "id") %>%
  distinct(id, .keep_all = TRUE) -> count_matrix_subset


# log2(RPKM)-transformation ----
data <- select(count_matrix_subset, -name, -pident, -qcovs)[-1] #restore before doing various data transformations
data <- log2(data+1)
data <- as.matrix(data) # convert to matrix
data <- t(scale(t(data), scale=F)) # Centering rows
data <- as.data.frame(data)

data$total <- log2(rowSums(select(count_matrix_subset, -name, -pident, -qcovs)[-1]))

total <-  log2(colSums(select(count_matrix_subset, -name, -pident, -qcovs)[-1]))

data <- data.frame(data, 
                   name = count_matrix_subset$name,
                   pident = count_matrix_subset$pident, 
                   qcovs = count_matrix_subset$qcovs)

rownames(data) <- data$name

point.col <- rep("grey", nrow(data))

#point.col[which(data$total >= round(mean(data$total)))]  <- 'red'
point.col[which(data$total == max(data$total))]  <- 'red'
highestp <- data[which(data$total == max(data$total)),'name']
# Viz parsed-expr ----
path <- "/Users/cigom/transcriptomics/miRNAs"
png(paste0(path,"/", "superheat.png"), units="px", width=2400, height=2600, res=400)

# superheat ----
superheat(dplyr::select(data, -total, -name, -pident, -qcovs),
          # add scatterplot next to the rows
          yr = data$total,
          yr.axis.name = paste0("log2(Expr)"), #, highestp, "\ncolored"),
          yr.plot.type = "scattersmooth",
          # change the color of the points
          yr.obs.col = point.col,
          # retain original order of rows/cols
          pretty.order.rows = TRUE,
          pretty.order.cols = TRUE,
          col.dendrogram = TRUE,
          # left labels
          left.label.size = 0.5,
          left.label.text.size = 2.7,
          left.label.text.alignment = 'right',
          #left.label.text.angle = 95,
          # bottom labels
          bottom.label.size = 0.5,
          bottom.label.text.size = 2.7,
          bottom.label.text.angle = 90,
          bottom.label.text.alignment = "right",
          row.title = " ",
          column.title = " ",
          # change the grid color
          grid.hline.col = "white",
          grid.vline.col = "white",
          left.label.col = "white",
          bottom.label.col = "white",
          heat.pal = brewer.pal(n = 9, name = "RdYlBu"))


dev.off()
# Selecto gonada data
gonada <- c("T0_R",  "T0_S1", "T1_R",  "T1_S2", "T2_R",  "T2_S3", "T3_R",  "T3_S4")
gonada_df <- select(data, gonada, total, name)


# Corr test----

library(corrplot)
library(dplyr)

# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}


gonada_cor <- t(select(gonada_df, -total, -name))
# matrix of the p-value of the correlation
p.mat <- cor.mtest(gonada_cor) #select(gonada_df, -total, -name))
# Correlation
corr <- cor(gonada_cor) #select(gonada_df, -total, -name))
# Viz corplot ----
corr[which(is.na(corr))] <- 0
p.mat[which(is.na(p.mat))] <- 0

corrplot(corr, method="color", type="upper", 
         order="hclust", 
         #col=brewer.pal(n=8, name="RdYlBu"),
         col=c("black", "white"),
         bg="lightblue",
         tl.srt=45, tl.col="black", tl.pos = "n",
         #p.mat = p.mat, sig.level = 0.01, insig = "blank",
         title = "\nCorplot of miRNA battery expressed in oyster samples (sig.level > 0.01 white color)")

str(gonada_df)
# BLASTN hairpin db vs genome----


blast_path = '/Users/cigom/transcriptomics/miRNAs/blastn'
blast_files <- list.files(blast_path, full.names = TRUE, pattern = 'hairpin_collapsed2dna_vs_GCA_002022765.4_C_virginica-3.0_genomic_blastn.outfmt6')

blast_plot <- do.call(rbind,lapply(blast_files, blast_load))

#GGally::ggpairs(blast_plot[,-c(1,2)])

# 
# hairpin_map_pos <- blast_plot[, c('qlen', 'qstart', 'qend')]
# genome_map_pos <- reshape2::melt(blast_plot[, c('slen', 'sstart', 'send')])

#
start <- blast_plot[, 'qstart']
end <- blast_plot[, 'qend']
df <- as.data.frame(cbind(start,end))

poss <- vector()

i=1
for(i in 1:nrow(df)){
  poss <- append(poss, c(df[i,1]:df[i,2]))
  i+1
}

dens <- density(poss)
plot(dens)


# Viz blastn ----
df <- dplyr::select(blast_plot, -qseqid, -sseqid)

p <- ggplot() +
  geom_density(data = reshape2::melt(df),
             mapping = aes(x = value,
                           y = ..count.. 
                           ), color = "grey") +
  labs(x = "Position", y = "Freq")

p + facet_wrap( ~ variable, scales = "free")


png(paste0(path,"/", "hairpin_collapsed2dna_vs_genomic_blastn.gapopen.png"), units="px", width=2700, height=2200, res=400)


p <- ggplot() +
     geom_point(data = blast_plot,
                mapping = aes(x = qcovs,
                              y = gapopen, color = gapopen <= 10), alpha = 0.7) +
     scale_colour_manual(values = c("black", "red")) +
     scale_x_continuous(name = "Coverage", limits =  c(25, 100)) +
     scale_y_continuous(name = 'gapopen') + geom_rug(data = blast_plot,
             mapping = aes(x = qcovs,
                           y = gapopen), alpha = 1/2, position = "jitter") +
  labs(
    title = 'Hairpin sequences aligned to the oyster genome',
    caption = 'hairpin_collapsed2dna_vs_GCA_002022765.4_C_virginica-3.0_genomic_blastn.outfmt6')

p
dev.off()
# 逆向きのhitは赤で表示する
# ref https://rstudio-pubs-static.s3.amazonaws.com/297114_2c17198cfa3c49ef87ad40e864cf6597.html
# library(devtools)
# install_github("trinker/pacman")

# blast distribution hits plot ----
ggblast <- function(blast_out, num_hit){
  # usage: ggblast(blastout=bnout, 　     # blast結果data.frame
  #                num_query=c(1,3,5),    # uniqueなqueryの番号、もしくはid
  #                num_hit=10             # top10
  #                 ) 
  pacman::p_load(RColorBrewer, zoo, dplyr, ggplot2)
  
  # Etiqueta de leyenda
  vlab <- c(paste0(zoo::rollapply(c(0,seq(70,100,5)), width=2, by=1, 
                                  function(x){paste(x, collapse = ":")}), "(-)"),
            zoo::rollapply(c(0,seq(70,100,5)), width=2, by=1, 
                           function(x){paste(x, collapse = ":")}))
  # identity - Codigo de color
  vcol <- c(brewer.pal(7, "Reds"), brewer.pal(7, "Blues"))
  if(missing(num_hit)){ num_hit <- sum(blast_out$qseqid %in% x) }
  # queryごとにplot -----
  gg_blast <-
    lapply(unique(blast_out$qseqid),
           function(x){
             # queryごとのデータフレーム
             qdat <-blast_out %>%
               # identityを離散値に変換
               mutate(pident=as.integer(cut(pident, breaks=c(0,seq(70,100,5)), 1:7))) %>%
               # hit position
               mutate(pident=ifelse((qstart-qend) > 0, pident, pident+7)) %>%
               # query selection
               filter(qseqid==x) %>%
               # top hitの数を指定
               filter(row_number() %in% 
                        1:ifelse(is.null(num_hit), sum(blast_out$qseqid %in% x), num_hit)) %>%
               # identity, hitをfactorにして水準を逆に
               mutate(pident=factor(pident, levels=rev(levels(factor(pident)))),
                      sseqid=factor(sseqid, levels=rev(unique(sseqid))))
             
             # queryごとのカラーコードとレジェンドラベル
             col_qdat <- vcol[as.integer(levels(factor(qdat$pident)))]
             lab_qdat <- vlab[as.integer(levels(factor(qdat$pident)))]
             
             # 
             ggplot(qdat, aes(x=qstart,xend=qend, y=sseqid, yend=sseqid, colour=pident)) +
               theme_bw() + labs(title=x) +
               geom_segment(size=3) +
               scale_color_manual(values = col_qdat, name='Key of\nAligment\nScore', labels = lab_qdat) +
               labs(title = paste0("Distribution of top ", num_hit, " blast hits")) + labs(x = 'Scafold Position (Genome)', y = 'Consensus hairpin sequence') +
               theme_classic(base_size = 12) 
           }
    )
  return(gg_blast)
}


blast_out <- blast_plot[blast_plot$gapopen <= 10,]
blast_out <- blast_plot[blast_plot$qend <= 100,]
blast_out <- blast_plot
#blast_out <- blast_plot
blast_out[,2] <- "Genome"
names(blast_out)[c(1,2)] <- c("sseqid","qseqid")

gg_bpout <- ggblast(blast_out = blast_out, num_hit = 100)

library(gridExtra)

#path <- "/Users/cigom/transcriptomics/miRNAs"
png(paste0(path,"/", "hairpin_collapsed2dna_vs_genomic_blastn..png"), units="px", width=2400, height=2600, res=400)

do.call(grid.arrange, gg_bpout) # gg_bpout[1:6]

dev.off()

blast_out$mmatch <- "> 5" 
# blast_out[which(blast_out$mismatch <= 20), 'mmatch'] <- '<=20'
# blast_out[which(blast_out$mismatch <= 10), 'mmatch'] <- '<=10'                     
blast_out[which(blast_out$mismatch <= 5), 'mmatch'] <- '<=5'  

blast_out$mmatch <- factor(blast_out$mmatch) #levels = c('1e-10', '1e-20', '1e-40', '1e-80', '1e-120', '1e-200', '1e-300'))

# 
ggplot(blast_out, aes(x=pident, y=length)) +
  geom_point(aes(color = gapopen, size = qcovs, shape = mmatch), alpha = 0.7) +
  scale_shape_manual(values = 14:15+length(levels(blast_out$mmatch))) + # if shape = name
  scale_color_distiller(palette = "YlOrBr") +
  theme_classic() +
  labs(x = 'Identidad del alineamiento (%)',
       y = 'Longitud del alineamiento',
       subtitle = 'Blast hits (blastn): Oyster genome vs precursor miRNA metazoa DataBase',
       caption = "")


# Load PFAM ----

#  tail -n +4 longest_orfs_PFAM.out | awk '{print $1,$2,$3,$6,$7,$8,$9,$24"_"$25, $4}' > longest_orfs_PFAM.headers.out

pfam_file = '/Users/cigom/transcriptomics/miRNAs/trinity_genes.fasta.transdecoder_dir/longest_orfs_PFAM.headers.out'

pfam_obj <- readr::read_delim(pfam_file, delim = " ", col_names = FALSE)

names(pfam_obj) <- c("target_name", "accession", "tlen", "qlen",   "E-value",  "score",  "bias", "description", "query_name")

pfam_obj$id <- sapply(strsplit(pfam_obj$query_name, "[::]"), `[`, 3)

# Argonaute select ----
# argonaute domains (pfam id):
# PF02171	Piwi	Piwi domain					
# PF02170	PAZ	PAZ domain


pfam_plot <-rbind(dplyr::filter(pfam_obj, grepl("PF02170",accession)),
            dplyr::filter(pfam_obj, grepl("PF02171",accession))) 


# Change point shapes, colors and sizes, shape=cyl,
ggplot(pfam_plot, aes(x=score, y=qlen)) +
  geom_point(aes(color = `E-value`, size = bias, shape = target_name), alpha = 0.7) +
  # scale_shape_manual(values = 1:length(levels(blastn_plot$name))) + # if shape = name
  scale_color_distiller(palette = "RdYlBu") +
  theme_classic() +
  labs(x = 'Puntaje',
       y = 'Longitud de secuencia',
       subtitle = 'per-domain hits',
       caption = "Oyster pseudo-genes againts a profile pfam domain DataBase (hmmscan)" )


# Parse to RPKM ----

pfam_plot %>% 
  group_by(id) %>%
  inner_join(count_matrix, by = "id") %>%
  distinct(id, .keep_all = TRUE) -> domain_exp

m <- reshape2::melt(select(domain_exp, names(count_matrix), target_name))

p <- ggplot(m, aes(x=variable, y= value, fill = target_name)) +
          geom_bar(stat="identity", color="black", position=position_dodge()) 
p + scale_fill_brewer(palette="Dark2") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Samples", y = "RPKM")


# Samples selection ----
# PCA ----
# select samples with highest variance in gene expression 
# to decipher whether RNAi is active

#data0 <- select(data, -total, -name) # if log2

data0 <- select(count_matrix_subset, -name, -pident, -qcovs)[-1] # if raw 

pca <- prcomp(t(data0), center = TRUE, scale. = TRUE)

## make a scree plot
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
# factoextra::fviz_pca_ind(pca) 
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

samples_data <- samples_data[match(rownames(pca$x), samples_data[,'sample_name']),]

# Or

if(identical(rownames(pca$x), samples_data[,'sample_name'])) { 
  pca.data <- data.frame(Sample=rownames(pca$x),
                         X=pca$x[,1],
                         Y=pca$x[,2],
                         Pollutant=samples_data[, 'Crude.oil..ug.L.'],
                         Tissue=samples_data[,'Tissue'])
 
  levels_hc <- c('0', '50', '100', '200')
  pca.data$Pollutant <- factor(pca.data$Pollutant, levels = levels_hc)
  pca.data$Tissue <- factor(pca.data$Tissue, levels = c('DG', 'GON'))
  
} else {
  pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
}



path <- "/Users/cigom/transcriptomics/miRNAs"
png(paste0(path,"/", "PCA.png"), units="px", width=2600, height=2200, res=400)

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
  #geom_text() +
  geom_point(aes(color = Tissue, shape = Pollutant), size = 5, alpha = 0.7) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  scale_shape_manual(values = 1:length(levels(pca.data$Pollutant)), name = "Pollutant (ug/L)") +
  scale_color_brewer(palette = "Set1") +
  ggtitle("samples explained by pca [miRNA battery subset]")

dev.off()

# and ANOVA ----

if(identical(colnames(data0), samples_data[,'sample_name'])) {
  sample_factoring <- factor(samples_data[,'Crude.oil..ug.L.'],
                              levels = c('0', '50', '100', '200'))
}

#

data0 = as.matrix(data0) # convert to matrix

anova_pvals = c()

for (j in 1:nrow(data0)) {
  feature_vals = data0[j,]
  data_for_anova = data.frame(y=feature_vals, group=factor(sample_factoring))
  fit = lm(y ~ group, data_for_anova)
  a = anova(fit)
  p = a$"Pr(>F)"[1]
  anova_pvals[j] = p
}


anova_ranking = order(anova_pvals)

fdr = p.adjust(anova_pvals, method = 'fdr')
anova_stats = data.frame(Pvals=anova_pvals, FDR=fdr)
rownames(anova_stats) = rownames(count_matrix_subset)
adj_data = cbind(data, anova_stats)
adj_data = adj_data[order(adj_data$Pvals),]
signif_indices = (adj_data$Pvals<=0.05) # instead of fdr use pvals

if (sum(signif_indices)==0) stop('No significant variable features identified. Stopping.');
adj_data = adj_data[signif_indices,]

rownames(adj_data) <- adj_data$name

point.col <- rep("grey", nrow(adj_data))
point.col[which(adj_data$FDR <= 0.05)] <- 'red'

png(paste0(path,"/", "superheat2.png"), units="px", width=2600, height=2800, res=400)

superheat(select(adj_data, -total, -name, -pident, -qcovs, -Pvals, -FDR),
          # add scatterplot next to the rows
          yr = adj_data$total,
          yr.axis.name = paste0("log2(rpkm)"), #\n[FDR <= 0.05 \ncolored]
          yr.plot.type = "scattersmooth",
          # change the color of the points
          yr.obs.col = point.col,
          # retain original order of rows/cols
          pretty.order.rows = TRUE,
          pretty.order.cols = TRUE,
          col.dendrogram = TRUE,
          # left labels
          left.label.size = 0.5,
          left.label.text.size = 2,
          left.label.text.alignment = 'right',
          #left.label.text.angle = 95,
          # bottom labels
          bottom.label.size = 0.5,
          bottom.label.text.size = 2,
          bottom.label.text.angle = 90,
          bottom.label.text.alignment = "right",
          row.title = " ",
          column.title = " ",
          # change the grid color
          grid.hline.col = "white",
          grid.vline.col = "white",
          left.label.col = "white",
          bottom.label.col = "white",
          heat.pal = brewer.pal(n = 9, name = "RdYlBu"))
# or
dev.off()

ggplot(adj_data, aes(x=pident, y=qcovs)) +
  geom_point(aes(color = Pvals , size = FDR), alpha = 0.7)

# 
data1 = adj_data[,colnames(adj_data) %in% colnames(data0)]
data1 = data1[1:min(nrow(data1),1000),] # restrict to 1000 with anova sig P-value ranking
data1 = data1[,colSums(data1)>0]

data1 = t(scale(t(data1), scale=F))
sample_dist = dist(t(data1), method='euclidean')
hc_samples = hclust(sample_dist, method='complete')

prin_comp_data = data1
pca = prcomp(t(prin_comp_data), center = FALSE, scale. = FALSE)
pc_pct_variance = (pca$sdev^2)/sum(pca$sdev^2)
barplot(pc_pct_variance, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

samples_data <- samples_data[match(rownames(pca$x), samples_data[,'sample_name']),]

library(factoextra)

# http://www.sthda.com/english/wiki/fviz-pca-quick-principal-component-analysis-data-visualization-r-software-and-data-mining

fviz_pca_ind(pca, label="none", habillage=samples_data$Crude.oil..ug.L.,
                  addEllipses=TRUE, ellipse.level=0.95) + scale_color_brewer(palette="Dark2") 

# and variance ----
# cite: https://davetang.org/muse/2012/04/14/variance-in-rna-seq-data/
dim(data0 <- select(count_matrix_subset, -name)[-1])
dim(data_subset <- data0[rowSums(data0) > 0,])

# prepare normalization library size due to:
colSums(data_subset)

#function for normalisation
norm_cpm <- function(x){ x * 1000000 / sum(x) }

#apply the normalisation
data_subset_cpm <- apply(data_subset, 2, norm_cpm)
head(data_subset_cpm)

# sanity check
colSums(data_subset_cpm)

data_var <- apply(data_subset_cpm, 1, var)
data_mean <- apply(data_subset_cpm, 1, mean)

par(mfrow=c(1,3))
plot(density(data_var))
plot(density(data_mean))
plot(log2(data_mean), log2(data_var), pch='.', main = paste0('Corr = ',cor(data_mean, data_var, method="spearman")))
#correlation of mean and variance; the higher the expression the higher than variance is

# Microtar ----
# miRNA prediction
# after install microtar v0.9.6
# ./microtar -t trinity_genes.fasta.transdecoder.cds -q mature.fa -f microtar_results.tsv


path_mirna <- '/Users/cigom/transcriptomics/miRNAs/mRNA_miRNA_prediction'
microtar_file <- list.files(path_mirna, pattern = ".tsv", full.names = TRUE)
microtar_obj <- read.csv(microtar_file, sep = '\t', header = FALSE)
names(microtar_obj) <- c("miRNA", "miLength", "gene", "mLength", "smatch", "energy", "norm_energy")

hist((microtar_obj$smatch))
plot(density((microtar_obj$energy)))
plot(microtar_obj$energy, microtar_obj$norm_energy)

summary(microtar_obj[microtar_obj$norm_energy > 1, c("smatch", "energy", "norm_energy")])

# mirDeep2 ----
# compare vs mirdeep2

# bowtie-build trinity_genes.fasta.transdecoder.pfam.cds trinity_genes.fasta.transdecoder.pfam.cds
mirdeep2_file <- list.files(path_mirna, pattern = "mature.headers_vs_GCA_002022765.4_C_virginica-3.0_genomic.fa.arf", full.names = TRUE)
mirdeep2_obj <- read.csv(mirdeep2_file, sep = '\t', header = FALSE, stringsAsFactors = FALSE)
names(mirdeep2_obj) <- c("readID_wo_whitespaces", "length_read", "start", "end", "read_sequence", "genomicID_wo_whitspaces", "length_genomic", "g_start", "g_end", "genomic_sequence", "strand", "mismatches", "editstring")
nrow(mirdeep2_obj)

# strand selection (+ / -) 
table(mirdeep2_obj$mismatches)

# N of mismatches
table(mirdeep2_obj[,12])

str(mirdeep2_obj)

# select 3pUTR ----
# convert bed 3p_UTR
UTR_3p_file <- list.files(path_mirna, pattern = "*bed", full.names = TRUE)
UTR_3p <- read.csv(UTR_3p_file, sep = '\t', header = FALSE)

UTR_3p$length <- UTR_3p[,5] - UTR_3p[,4]

UTR_3p_out <- UTR_3p[UTR_3p$length >= 100,]

hist(UTR_3p_out$length, breaks = 1000, 
                        main = paste0("the minLenght: ", min(UTR_3p_out$length), "\nthe maxLenght: ", max(UTR_3p_out$length)))

write.table(UTR_3p_out[, -c(ncol(UTR_3p_out))], file =  paste0(path_mirna, "/", "three_prime_UTR.bed"), 
                            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')

# miRanda ----
path_miranda <- c('/Users/cigom/transcriptomics/miRNAs/mRNA_miRNA_prediction/miranda')
miRanda_file <- list.files(path_miranda, pattern = "mature_collapsed2dna_vs_trinity_genes.3p_UTR.R.txt", full.names = TRUE)
miRanda_obj <- read.csv(miRanda_file, sep = '\t', header = FALSE, stringsAsFactors = FALSE)

names <- c("Seq1","Seq2","TotScore","TotEnergy","MaxScore","MaxEnergy","Strand","Len1","Len2","Positions")

names(miRanda_obj) <- names
str(miRanda_obj)

breaks <- c(-60, -40, -35, -30, 25)
breaks <- c(round(summary(miRanda_obj$MaxEnergy))) - 10
names(breaks) <- breaks

p <- ggplot() +
  stat_summary_2d(data = miRanda_obj,
                     mapping = aes(x = Len2,
                                   y = MaxScore,
                                   z = MaxEnergy),
                              fun = mean) +
      # scale_fill_gradient(name = "Folding Energy\n(kcal/mol)",
      #                 low = "black",
      #                 high = "#BFBCFF",
      #                 breaks = breaks, 
      #                 labels = format(breaks)) +
  labs(x='trinity_genes 3p-UTR length')

# WGCNA ----
# # #
# install.packages("BiocManager")
# BiocManager::install("WGCNA")

library("WGCNA")

# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html
# https://pdfs.semanticscholar.org/b756/a0c832c4ec3fa6c263026a9b025325e2abf8.pdf
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

datExpr0 = as.data.frame(t(count_matrix[-ncol(count_matrix)]))
names(datExpr0) = count_matrix$id
rownames(datExpr0) = names(count_matrix[, -c(1)])

gsg = goodSamplesGenes(datExpr0, verbose = 3)

# Excluding 2264 genes from the calculation due to too many missing samples or zero variance
# If the last statement returns TRUE, all genes have passed the cuts. 
# If not, we remove the offending genes and samples from the data:

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))

plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


# check outliers to remove:
# h <- max(sampleTree$height)
# sampleTree[sampleTree$height == max(sampleTree$height), 'labels']
# then, remove:

# Plot a line to show the cut
abline(h = 3e+05, col = "red");

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 3e+05, minSize = 0)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# ^ The variable datExpr now contains the expression data ready for network analysis.


datExpr_subset <- datExpr[names(datExpr) %in% count_matrix_subset$id]
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datExpr_subset, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = count_matrix_subset$name,
                    main = "Sample dendrogram and trait heatmap")


# One-step network construction and module detection

# enableWGCNAThreads()

net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)

# We now return to the network analysis. 
# To see how many modules were identified and what the module sizes are, one can use table(net$colors). Its output is
table(net$colors)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, mergedColors)$eigengenes

# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5)

par(cex = 0.9)

plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), 
                      marHeatmap = c(3,4,1,2), 
                      cex.lab = 0.8, xLabelsAngle = 90)

# auto-curated sequence annotation with ARBitrator for miRNA battery
# superiority ----
# cd02040 (nifH) = determined (curated) target
# We define a sequence’s superiority as log 10 ( E -value of best non-nhit) −log 10 ( E -value of cd02040 hit).
# Thus, candidates are accepted if their superiority is ≥1. 
# The output of the sensitivity phase is a subset of the protein GIs generated by the sensitivity phase, representing sequences that ARBitrator classifies as nifH 
# superiority measures the degree to which a candidate is more similar to nifH than to genes coding for other known proteins
# https://academic.oup.com/bioinformatics/article/30/20/2883/2422235
# Quality ----
# For all candidate sequences (blast hits) retrieved, 
# we define ‘quality’ as the negative log 10 of the E -value of the hit. 
# If a subject is hit by multiple queries from the representative set, 
# then quality is defined as the negative log 10 of the smallest E -value across all hits. 
# The sensitivity phase accepts sequences with quality ≥2 (i.e. all E -values are ≤0.01)