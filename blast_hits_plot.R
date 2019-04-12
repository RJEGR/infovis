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

blast_path = '/Users/cigom/transcriptomics/miRNAs/blastp'
blast_files <- list.files(blast_path, full.names = TRUE, pattern = 'trinity_genes.fasta.transdecoder_vs_peerj-04-1763-s005_blastp.outfmt6')


blast_load <- function(x) { 
  blast_obj <- read.table(x, stringsAsFactors = FALSE)
  blast_names <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs", "qcovhsp", "qlen", "slen")
  names(blast_obj) <- blast_names
  #blast_obj <- blast_obj[, c("pident", "mismatch","gapopen", "qcovs", "length")]
  return(blast_obj)
}

blast_plot <- do.call(rbind,lapply(blast_files, blast_load))

#blastn_plot$name <- factor(blastn_plot$name, levels = c('1e-10', '1e-20', '1e-40', '1e-80', '1e-120', '1e-200', '1e-300'))
n_sseqid <- 10
ssqid_text <- names(table(blast_plot$sseqid)[table(blast_plot$sseqid) > n_sseqid])
plot_text <- blast_plot[blast_plot$pident >= 90 & blast_plot$qcovs >= 90, 'sseqid']
library(ggplot2)

# Change point shapes, colors and sizes, shape=cyl,
ggplot(blast_plot, aes(x=pident, y=length)) +
  geom_point(aes(color = qcovs, size = mismatch), alpha = 0.7) +
  # scale_shape_manual(values = 1:length(levels(blastn_plot$name))) + # if shape = name
  scale_color_distiller(palette = "RdYlBu") +
  theme_classic() +
  labs(x = 'Identidad del alineamiento (%)',
       y = 'Longitud del alineamiento',
       subtitle = 'Blast hits (blastp): Oyster pseudo-genes vs miRNA machinery DataBase',
       #caption = paste0('miRNA machinery with Identity >= 90 & Coverage >= 90 is labeled')) +
       caption = paste0('miRNA machinery with Frequency > ', n_sseqid, " is labeled")) +
  #geom_text(data = subset(blast_plot, pident >= 90 & qcovs >= 90 ), aes(label=sseqid), hjust = 1, vjust = 1, check_overlap =TRUE, size = 2.3)
  geom_text(data = subset(blast_plot, sseqid==ssqid_text), aes(label=sseqid), hjust = 0, vjust = 1, check_overlap =TRUE, size = 2.3)

  #xlim(c(min(blastn_plot$qlen), max(blastn_plot$qlen)))
  #facet_wrap( ~ name, scales = 'free')


# PFAM 
# awk '{print $1}' longest_orfs_PFAM.out | sort | uniq -c | sort -k1,1 | awk '{print $2, $1}' > longest_orfs_PFAM.out.freq

pfam_file = '/Users/cigom/transcriptomics/miRNAs/trinity_genes.fasta.transdecoder_dir/longest_orfs_PFAM.out.freq'
pfam_obj <- read.table(pfam_file, stringsAsFactors = FALSE)

plot(density(pfam_obj$V2))


