
path <- '/Users/cigom/transcriptomics/oktopus_full_assembly/all_deg_ricardo'

files <- list.files(path = path, pattern = 'annotation.txt', full.names = T)

# fileName <- files[1]


x <- lapply(files, function(x) {
  
  group <- paste(sapply(strsplit(basename(x), "[.]"), `[`, 1),collapse = '_')

  y <- read.csv(x, sep = '\t')
  y <- data.frame(y, group = group)
  return(y)})

x <- do.call(rbind, x) %>% as_tibble()

query.genes <- x %>% distinct(ID) %>% pull()

str(query.genes)


library(tidyverse)
library(topGO)

runtopGO <- function(topGOdata, topNodes = 20, conservative = TRUE) {
  
  RFisher <- runTest(topGOdata, 
                     algorithm = "classic", 
                     statistic = "fisher")
  
  # To make this test conservative. Next we will test the enrichment using the Kolmogorov-Smirnov test. We will use the both the classic and the elim method.
  
  if(conservative) 
  {
    RKS <- runTest(topGOdata, algorithm = "classic", 
                   statistic = "ks")
    
    RKS.elim <- runTest(topGOdata, algorithm = "elim", 
                        statistic = "ks")
    
    
    
    
    allRes <- GenTable(topGOdata, 
                       classicFisher = RFisher,
                       classicKS = RKS, 
                       elimKS = RKS.elim,
                       orderBy = "elimKS", 
                       ranksOf = "classicFisher", 
                       topNodes = topNodes) 
  } else {
    RKS <- runTest(topGOdata, algorithm = "classic", 
                   statistic = "ks")
    
    test.stat <- new("weightCount",
                     testStatistic = GOFisherTest,
                     name = "Fisher test", sigRatio = "ratio")
    
    weights <- getSigGroups(topGOdata, test.stat)
    
    allRes <- GenTable(topGOdata,
                       classic = RFisher,
                       KS = RKS,
                       weight = weights,
                       orderBy = "weight",
                       ranksOf = "classic",
                       topNodes = topNodes)
    
    # allRes <- GenTable(object = topGOdata, 
    #                    elimFisher = RFisher,
    #                    topNodes = topNodes)
  }
  
  return(allRes)
}

description <- "complete topGO enrichment using split_annot"

path <- "~/transcriptomics/oktopus_full_assembly/"

go_file <- paste0(path, '/Trinotate.xls.gene_ontology')

MAP <- topGO::readMappings(go_file)

path_out <- paste0(path, 'ANNOT_out')

# load(paste0(path, 'splited_annot.RData'))
# split_annot <- readRDS(paste0(path, 'biological_process.rds'))

 #...
# ya esta hecho el enriquecimiento de genes /Users/cigom/transcriptomics/oktopus_full_assembly/Enriquecimiento_GO
#  pero no sirve del todo, 

path <- '/Users/cigom/transcriptomics/oktopus_full_assembly/Enriquecimiento_GO/IDs_dentro_cat_GO/'

files <- list.files(path = path, pattern = 'txt', full.names = T)

x <- lapply(files, function(x) {
  
  group <- paste(sapply(strsplit(basename(x), "[.]"), `[`, 1),collapse = '_')
  
  y <- data.table::fread(x, header = F)
  y <- data.frame(y, group = group)
  return(y)})

x <- do.call(rbind, x) %>% as_tibble()
