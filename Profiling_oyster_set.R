#!/usr/bin/env Rscript

# usar clusterprofile sobre los genes expresados diff y la informacion de la anotacion Trinotate.xls

# ========
# Reading Data
# ========

rm(list=ls())
options(stringsAsFactors = FALSE)
# ========

set.seed(20200717)

.LoadToEnvironment <- function(RData, env = new.env()){
  load(RData, env)
  return(env)
}

print('Reading files.')

path <- "~/Documents/paper_oyster_2020/"

Rdata_file <- list.files(path = path,
                         pattern = 'Trinotate.xls.RData',
                         full.names = T)

count_file <- list.files(path = path,
                         pattern = 'RSEM.isoform.counts.matrix', 
                         full.names = T)

annot.env <- .LoadToEnvironment(Rdata_file)


y <- read.table(count_file, header=T, com='', 
                row.names=1, check.names=F, sep='\t', 
                stringsAsFactors = FALSE)

DE_file <- list.files(path = paste0(path, '/Up_down_list/'),
                      pattern = 'UP.txt', 
                      full.names = T)

read.group.table <- function(file) {
  read.table(file) %>%
    data.frame(., group = basename(file))
}

dim(de_df <- do.call(rbind, lapply(DE_file, read.group.table)))

de_df %>%
  mutate(group = sapply(strsplit(group, "[.]"), "[", 1)) %>%
  dplyr::rename(ids = V1) -> de_df
  

dim(data <- round(y) )
dim(data <- data[rowSums(edgeR::cpm(data) > 1) >= 2,])

data <- log2(data+1)
data <- as.matrix(data) # convert to matrix
data <- t(scale(t(data), scale=F)) # Centering rows
data <- as.data.frame(data)

# ========================
# Defining variables (some)
# ========================

go <- annot.env$go
blastx <- annot.env$blastx
egg <- annot.env$egg

# blastp <- annot.env$blastp
# pfam <- annot.env$pfam
# x <- annot.env$x

outpath <- paste0(path, '/profiling')
system(paste0('mkdir -p ', outpath))

# ===============
# Loading package
# ================

.cran_packages <- c('tidyverse', 'purrr', 
                    'tibble', 'reshape2', 'ggplot2')

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
   install.packages(.cran_packages[!.inst], dep=TRUE, repos='http://cran.us.r-project.org')
}

# Load packages into session, and print package version

sapply(c(.cran_packages), require, character.only = TRUE)

# prepare data-base reference ----

library(biomaRt)

# http://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html

# biomaRt::useEnsembl()
# biomaRt::listEnsemblArchives()
# biomaRt::listEnsembl()

ensembl_metazoa <- useEnsembl(biomart = "metazoa_mart",
                              host = "metazoa.ensembl.org")

# useMart(biomart = "metazoa_mart",
#         host = "metazoa.ensembl.org")

datasets <- listDatasets(ensembl_metazoa)

# head(datasets)

datasets %>%
  filter(str_detect(description, 'gigas'))

# searchDatasets(mart = ensembl_metazoa, pattern = "cgigas")

ensembl_cgigas <- useEnsembl(biomart = "metazoa_mart",
                             host = "metazoa.ensembl.org",
                             dataset = 'cgigas_eg_gene')

# =================
# Filtering dataset
# =================

table(go$ontology)

# The first two columns are: “gene_id” and “transcript_id”, representing predicted genes and their corresponding transcripts, respectively. The columns “Top_BLASTX_hit” and “Top_BLASTP_hit” show the top BlastX and BlastP hit results of homology searches against the NCBI database. BlastX is one of the latest additions to the Trinotate annotation pipeline and compares all six open reading frames (ORF) of the query sequences against the protein database. The RNAMMER column shows information about predicted ribosomal RNA genes discovered in the transcriptome assembly that were predicted by hidden Markov models (HMM). The prot_id, prot_coords and prot_seq columns provide the ID, location and translation of the longest ORFs, respectively. The Pfam column represents the HMMER/PFAM protein domain identification search results. HMMER is used to search databases for homologs of proteins, employing hidden Markov models. The SignalP column shows the presence and location of predicted signal peptides. Similarly, the TmHMM column presents the predicted transmembrane regions. The eggNOG (Evolutionary genealogy of genes: Non-supervised Orthologous Groups) column has the search result of the database of orthologous groups of genes, which are further annotated with functional description lines. Lastly, the gene_ontology column shows the relationship of our data to the Gene Ontology (GO) terms that aim to unify the representation of genes and gene products across all species ----



# continue ----

T1_100 <- 'T1_100'
T1_200 <- 'T1_200'
T3_100 <- 'T3_100'
T3_200 <- 'T3_200'

# GOType <- 'biological_process' # 5773 DiffEx
# GOType <- 'molecular_function' # 2958 DiffEx
# GOType <- 'cellular_component' # 3534 DiffEx


getList <- function(go, de_df, f_pattern, GOtype = 'MF') {
  
  term_type <- c('MF' = 'molecular_function', 
                 'CC' = 'cellular_component', 
                 'BP' = 'biological_process')
  
  GOtype <- toupper(GOtype)
  
  which_term <- which(GOtype %in% names(term_type))
  
  GOtype <- term_type[which_term]

  
  go %>%
    ungroup() %>%
    inner_join(de_df, by = c("transcript"='ids')) %>%
    filter_all(any_vars(str_detect(., GOtype))) %>%
    ungroup() %>%
    filter(str_detect(group, paste0(f_pattern, '_UP'))) %>%
    mutate(ontology = factor(ontology)) %>%
    group_by(transcript)
}



getList(go,de_df, T1_100, GOtype='MF') # 448/727 in GOType (MF)
getList(go,de_df, T1_200, GOtype='MF') # 224/407 in GOType (MF)
getList(go,de_df, T3_100, GOtype='MF') # 22/53 in GOType (MF)
getList(go,de_df, T3_200, GOtype='MF') # 236/380 in GOType (MF)


# table(de_df$group)


# Select terms from the universe 


# filters = listFilters(ensembl_cgigas)
# filters %>% filter(str_detect(name, 'go'))
# go_parent_term

# attributes = listAttributes(ensembl_cgigas)

# attributes %>% filter(str_detect(description, 'term'))

# enrichement by GO/entrezid

getEntrez <- function(go, de_df, pattern = NULL, GOtype='MF', mart) {
  
  library(biomaRt)
  
  if(is.null(pattern)) {
    geneList <- getList(go, de_df, '', GOtype)
  } else
    geneList <- getList(go, de_df, pattern, GOtype)
  
  
  
  GO_universe <- geneList %>% select_at(vars(go)) %>% pull(.)
  
  goids = getBM(attributes = c('entrezgene_id', 
                               'go_id',
                               'name_1006'),
                filters = 'go', 
                values = GO_universe, 
                mart = mart)
  
  
  entrezgene_id <- goids %>%
    filter(!is.na(entrezgene_id)) %>%
    mutate(entrezgene_id = as.character(entrezgene_id)) %>%
    arrange(desc(entrezgene_id))
  
  return(entrezgene_id)
}

# Example
getEntrez(go, de_df, GOtype='MF', mart = ensembl_cgigas,
          pattern = NULL) -> mart_data



enrichment <- function(go, de_df, GOtype='MF', pattern, mart = NULL) {
  
  library(clusterProfiler)
  library(enrichplot)
  
  geneList <- getList(go,de_df, pattern, GOtype)
  
  if(is.null(mart)) {
    
    ensembl_mart <- biomaRt::useEnsembl(biomart = "metazoa_mart",
                                 host = "metazoa.ensembl.org",
                                 dataset = 'cgigas_eg_gene')
    
    entrezgene_id <- getEntrez(go, de_df, pattern, GOtype, 
                               ensembl_mart)

  } else
    entrezgene_id <- mart
  
  
  entrezgene_id <- entrezgene_id %>%
    distinct() %>%
    arrange(desc(entrezgene_id))
  
  entrezgene_id %>%
    pull(entrezgene_id) -> gene_universe

  geneList %>%
    inner_join(entrezgene_id,
               by = c('go' = 'go_id')) %>% 
    distinct(entrezgene_id, go, transcript) %>%
    ungroup() %>%
    arrange(desc(entrezgene_id)) %>%
    pull(entrezgene_id) ->  gene
  
  TERM2GENE <- entrezgene_id[, c(2,1)]
  TERM2NAME <- entrezgene_id[, c(2,3)]
  
  ego <- enricher(gene =  gene,
                  universe = gene_universe,
                  pAdjustMethod = "fdr",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,
                  minGSSize = 10,
                  TERM2GENE = TERM2GENE,
                  TERM2NAME = TERM2NAME)
  
  # visualization ----
  
  plot_gsea <- function(enrichResult) {
    
    enrichResult %>%
      as.data.frame() %>%
      separate(GeneRatio, c('aRatio','bRatio')) %>%
      mutate(GeneRatio = as.numeric(aRatio) / as.numeric(bRatio)) %>%
      arrange(desc(GeneRatio)) %>%
      mutate(Description = factor(Description, levels = unique(Description))) %>%
      ggplot(aes(fill = p.adjust, color = p.adjust)) +
      # geom_bar(aes(y = Count, x = reorder(Description, Count)), width=0.7, stat = "identity") +
      geom_point(aes(size = Count, y = GeneRatio, 
                     x = reorder(Description, Count))) +
      ggsci::scale_fill_gsea() +
      ggsci::scale_color_gsea() +
      coord_flip() +
      labs(x = 'Description') +
      scale_size(range=c(3, 8)) +
      theme_bw() + 
      ggthemes::theme_tufte(base_family='GillSans') +
      theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                       size = 12, hjust = 1),
            axis.text.y = element_text(size = 12)) +
      guides(fill = guide_colorbar(barheight = unit(4, "in"),
                                   ticks.colour = "black", 
                                   frame.colour = "black",
                                   label.theme = element_text(size = 12)))
  }
  
  # barplot ----
  barplot_ego <- barplot(ego) +
    theme_bw() + 
    labs(x = 'Description', y = 'Count') +
    ggthemes::theme_tufte(base_family='GillSans') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 12, hjust = 1),
          axis.text.y = element_text(size = 12)) +
    ggsci::scale_fill_gsea() +
    ggsci::scale_color_gsea() +
    guides(fill = guide_colorbar(barheight = unit(4, "in"),
                                 ticks.colour = "black", 
                                 frame.colour = "black",
                                 label.theme = element_text(size = 12)))
  
  # emaplot ----
  emapplot_ego <- emapplot(ego) +
    theme_bw() + 
    ggthemes::theme_tufte(base_family='GillSans') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 12, hjust = 1),
          axis.text.y = element_text(size = 12)) +
    ggsci::scale_fill_gsea() +
    ggsci::scale_color_gsea() +
    guides(fill = guide_colorbar(barheight = unit(4, "in"),
                                 ticks.colour = "black", 
                                 frame.colour = "black",
                                 label.theme = element_text(size = 12)))
  # save ----
  
  ggsave(plot_gsea(ego),
         filename = paste0(pattern,"_plot_gsea.png"),
         path = outpath,
         width = 14, height = 12)
  
  ggsave(barplot_ego,
         filename = paste0(pattern,"_barplot.png"),
         path = outpath,
         width = 14, height = 7)
  
  ggsave(emapplot_ego,
         filename = paste0(pattern,"_emapplot.png"),
         path = outpath,
         width = 14, height = 7)
  
  return(ego)
}

# CC, MF, BP,

term <- 'BP'

T1_100_ego <- enrichment(go, de_df, GOtype = term,
                         pattern = T1_100, 
                         mart = mart_data)

T1_200_ego <- enrichment(go, de_df, GOtype = term,
                         pattern = T1_200,
                         mart = mart_data)

T3_100_ego <- enrichment(go, de_df, GOtype= term,
                         pattern = T3_100, 
                         mart = mart_data)

T3_200_ego <- enrichment(go, de_df, GOtype= term,
                         pattern = T3_200,
                         mart = mart_data)


# Kegg -----

keggEnrich <- function(go, de_df, pattern, GOtype, mart = NULL) {
  
  geneList <- getList(go, de_df, pattern, GOtype)
  
  if(is.null(mart)) {
    
    ensembl_mart <- biomaRt::useEnsembl(biomart = "metazoa_mart",
                                        host = "metazoa.ensembl.org",
                                        dataset = 'cgigas_eg_gene')
    
    entrezgene_id <- getEntrez(go, de_df, pattern, GOtype, 
                               ensembl_mart)
    
  } else
    entrezgene_id <- mart
  
  
  entrezgene_id <- entrezgene_id %>%
    distinct() %>%
    arrange(desc(entrezgene_id))
  
  # entrezgene_id %>%
  #   pull(entrezgene_id) -> gene_universe
  
  geneList %>%
    inner_join(entrezgene_id,
               by = c('go' = 'go_id')) %>% 
    distinct(entrezgene_id, go, transcript) %>%
    ungroup() %>%
    arrange(desc(entrezgene_id)) %>%
    pull(entrezgene_id) %>%
    bitr_kegg(., fromType = "kegg",
              toType = "Path", organism = "crg") %>%
    pull(kegg) %>%
    enrichKEGG(.,organism     = 'crg',
               pvalueCutoff = 0.05) -> keggE
  
  
  plot_gsea <- function(enrichResult) {
    
    enrichResult %>%
      as.data.frame() %>%
      separate(GeneRatio, c('aRatio','bRatio')) %>%
      mutate(GeneRatio = as.numeric(aRatio) / as.numeric(bRatio)) %>%
      arrange(desc(GeneRatio)) %>%
      mutate(Description = factor(Description, levels = unique(Description))) %>%
      ggplot(aes(fill = p.adjust, color = p.adjust)) +
      # geom_bar(aes(y = Count, x = reorder(Description, Count)), width=0.7, stat = "identity") +
      geom_point(aes(size = Count, y = GeneRatio, 
                     x = reorder(Description, Count))) +
      ggsci::scale_fill_gsea() +
      ggsci::scale_color_gsea() +
      coord_flip() +
      labs(x = 'Description') +
      scale_size(range=c(3, 8)) +
      theme_bw() + 
      ggthemes::theme_tufte(base_family='GillSans') +
      theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                       size = 12, hjust = 1),
            axis.text.y = element_text(size = 12)) +
      guides(fill = guide_colorbar(barheight = unit(4, "in"),
                                   ticks.colour = "black", 
                                   frame.colour = "black",
                                   label.theme = element_text(size = 12)))
  }
  
  # barplot ----
  barplot <- barplot(keggE) +
    theme_bw() + 
    labs(x = 'Description', y = 'Count') +
    ggthemes::theme_tufte(base_family='GillSans') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 12, hjust = 1),
          axis.text.y = element_text(size = 12)) +
    ggsci::scale_fill_gsea() +
    ggsci::scale_color_gsea() +
    guides(fill = guide_colorbar(barheight = unit(4, "in"),
                                 ticks.colour = "black", 
                                 frame.colour = "black",
                                 label.theme = element_text(size = 12)))
  
  # emaplot ----
  emapplot <- emapplot(keggE) +
    theme_bw() + 
    ggthemes::theme_tufte(base_family='GillSans') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 12, hjust = 1),
          axis.text.y = element_text(size = 12)) +
    ggsci::scale_fill_gsea() +
    ggsci::scale_color_gsea() +
    guides(fill = guide_colorbar(barheight = unit(4, "in"),
                                 ticks.colour = "black", 
                                 frame.colour = "black",
                                 label.theme = element_text(size = 12)))
  # save ----
  
  ggsave(plot_gsea(keggE),
         filename = paste0(pattern,"_barplot_kegg.png"),
         path = outpath,
         width = 14, height = 12)
  
  ggsave(barplot,
         filename = paste0(pattern,"_barplot_kegg.png"),
         path = outpath,
         width = 14, height = 7)
  
  ggsave(emapplot,
         filename = paste0(pattern,"_emapplot_kegg.png"),
         path = outpath,
         width = 14, height = 7)
  
}

# go, de_df, pattern, GOtype, mart = NULL
T1_100_kegg <- keggEnrich(go, de_df, T1_100,
                          GOtype = term,
                          mart = mart_data)

T1_200_kegg <- keggEnrich(go, de_df, T1_200,
                          GOtype = term,
                          mart = mart_data)

T3_100_kegg <- keggEnrich(go, de_df, T3_100,
                          GOtype = term,
                          mart = mart_data)

T3_200_kegg <- keggEnrich(go, de_df, T3_200,
                          GOtype = term,
                          mart = mart_data)


# done here!



