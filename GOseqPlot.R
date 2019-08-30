
rm(list=ls()); 


library(goseq)
library(GO.db)
library(qvalue)

dir <- '/Users/cigom/transcriptomics/oyster_gonada'
setwd(dir)

# capture list of genes for functional enrichment testing
factor_labeling = read.table("factor_labeling.txt", row.names=2, header=F)
colnames(factor_labeling) = c('type')
factor_list = unique(factor_labeling[,1])
DE_genes = rownames(factor_labeling)


# get gene lengths
gene_lengths = read.table("Trinity.fasta.seq_lens", header=T, row.names=1, com='')
gene_lengths = as.matrix(gene_lengths[,1,drop=F])


# get background gene list
background = read.table("background.txt", header=FALSE, row.names=1)
background.gene_ids = rownames(background)
background.gene_ids = unique(c(background.gene_ids, DE_genes))
sample_set_gene_ids = background.gene_ids


# parse GO assignments
GO_info = read.table("go_annotations.txt", header=F, row.names=1,stringsAsFactors=F)

GO_info_listed = apply(GO_info, 1, function(x) unlist(strsplit(x,',')))
names(GO_info_listed) = rownames(GO_info)
get_GO_term_descr =  function(x) {
    d = 'none';
    go_info = GOTERM[[x]];
    if (length(go_info) >0) { d = paste(Ontology(go_info), Term(go_info), sep=' ');}
    return(d);
}


#organize go_id -> list of genes
# time demand
GO_to_gene_list = list()

for (gene_id in intersect(names(GO_info_listed), sample_set_gene_ids)) {
    go_list = GO_info_listed[[gene_id]]
    for (go_id in go_list) {
        GO_to_gene_list[[go_id]] = c(GO_to_gene_list[[go_id]], gene_id)
    }
}


# GO-Seq protocol: build  probability weighting function (pwf) based on ALL DE features
# The PWF quantifies how the probability of a gene selected as DE changes as a function of its transcript length.
sample_set_gene_lengths = gene_lengths[sample_set_gene_ids,]
GO_info_listed = GO_info_listed[ names(GO_info_listed) %in% sample_set_gene_ids ]
cat_genes_vec = as.integer(sample_set_gene_ids %in% rownames(factor_labeling))
pwf=nullp(cat_genes_vec, bias.data=sample_set_gene_lengths)
rownames(pwf) = sample_set_gene_ids

plotPWF(pwf, binsize=1000)

# perform functional enrichment testing for each category.
for (feature_cat in factor_list) {
    message('Processing category: ', feature_cat)
    gene_ids_in_feature_cat = rownames(factor_labeling)[factor_labeling$type == feature_cat]
    cat_genes_vec = as.integer(sample_set_gene_ids %in% gene_ids_in_feature_cat)
    pwf$DEgenes = cat_genes_vec
    res = goseq(pwf,gene2cat=GO_info_listed, use_genes_without_cat=TRUE)
    
    ## over-represented categories:
     pvals = res$over_represented_pvalue
     pvals[pvals > 1 - 1e-10] = 1 - 1e-10
     q = qvalue(pvals)
     res$over_represented_FDR = q$qvalues
    go_enrich_filename = paste(feature_cat,'.GOseq.enriched', sep='')
    result_table = res[res$over_represented_pvalue<=0.05,]
    descr = unlist(lapply(result_table$category, get_GO_term_descr))
    result_table$go_term = descr;
    result_table$gene_ids = do.call(rbind, lapply(result_table$category, function(x) { 
            gene_list = GO_to_gene_list[[x]]
            gene_list = gene_list[gene_list %in% gene_ids_in_feature_cat]
            paste(gene_list, collapse=', ');
     }) )
    write.table(result_table[order(result_table$over_represented_pvalue),], file=go_enrich_filename, sep='	', quote=F, row.names=F)
    ## under-represented categories:
     pvals = res$under_represented_pvalue
     pvals[pvals>1-1e-10] = 1 - 1e-10
     q = qvalue(pvals)
     res$under_represented_FDR = q$qvalues
    go_depleted_filename = paste(feature_cat,'.GOseq.depleted', sep='')
    result_table = res[res$under_represented_pvalue<=0.05,]
    descr = unlist(lapply(result_table$category, get_GO_term_descr))
    result_table$go_term = descr;
    write.table(result_table[order(result_table$under_represented_pvalue),], file=go_depleted_filename, sep='	', quote=F, row.names=F)

}


quit(save = 'no')
# Manually plot
# Plot top 10 of GO enrichment

goseq_ <- function(pwf, factor_labeling, sample_factor, background) {
  gene_ids_in_feature_cat = rownames(factor_labeling)[factor_labeling$type == sample_factor]
  cat_genes_vec = as.integer(background %in% gene_ids_in_feature_cat)
  pwf$DEgenes = cat_genes_vec
  res = goseq(pwf,gene2cat=GO_info_listed, use_genes_without_cat=TRUE)
  return(res)
}

T0_res <- goseq_(pwf, factor_labeling, sample_factor = 'T0', background = sample_set_gene_ids)
T3_res <- goseq_(pwf, factor_labeling, sample_factor = 'T3', background = sample_set_gene_ids)

T0_res$sample <- 'T0'
T3_res$sample <- 'T3'

res <- rbind(T0_res, T3_res)

# https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
library(dplyr)
library(ggplot2)


## visualice:

res <- read.csv('T3.GOseq.enriched', header = TRUE, row.names = 1, sep = '\t')


res %>% 
  top_n(25, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  na.omit() %>%
  ggplot(aes(x=hitsPerc, y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) + geom_point() + expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count", 
       title = paste('GOseq enriched', sep='')) + theme_bw() +
  facet_wrap(~ sample)

#GOTERM[[res$category[1]]]

# The excess of low P-values shows that there are many GO 
# categories that contain a set of significantly long or short genes.


# Then, use semantic analysis.



