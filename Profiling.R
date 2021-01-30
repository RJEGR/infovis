#!/usr/bin/env Rscript

# ========
# Reading Data
# ========

rm(list=ls())

# ========

args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop("!!!\n A MATRIX OF DGE AND RDATA FROM ANNOTATION STEP IS NEEDED. 
         PLEASE, INPUT BOTH FILES IN THE SYNTAXIS AS FOLLOW EXAMPLE:\n
         Rscript --vanilla Profiling.R Trinotate.xls.RData diffExpr.matrix .\n", call.=FALSE)
} else {
  .LoadToEnvironment <- function(RData, env = new.env()){
  load(RData, env)
  return(env)  # # https://www.r-bloggers.com/safe-loading-of-rdata-files-2/
}
  print('Reading files.')
  annot.env <- .LoadToEnvironment(args[1])
  file = args[2]
  y <- read.table(file, header=T, com='', row.names=1, check.names=F, sep='\t', stringsAsFactors = FALSE)
  data <- y # restore before doing various data transformations
  data <- log2(data+1)
  data <- as.matrix(data) # convert to matrix
  data <- t(scale(t(data), scale=F)) # Centering rows
  data <- as.data.frame(data)
  data$transcript <- rownames(y)

  
}
# ========================
# Defining variables (some)
# ========================

go <- annot.env$go
blastx <- annot.env$blastx
blastp <- annot.env$blastp
pfam <- annot.env$pfam
x <- annot.env$x

outpath <- getwd()

# ===============
# Loading package
# ================

.cran_packages <- c('dplyr', 'purrr', 'tibble', 'reshape2', 'ggplot2')

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
   install.packages(.cran_packages[!.inst], dep=TRUE, repos='http://cran.us.r-project.org')
}

# Load packages into session, and print package version
sapply(c(.cran_packages), require, character.only = TRUE)

# =================
# Filtering dataset
# =================

#data = data.frame(transcript=sample(rownames(data), 10, replace = TRUE))

blastx %>%
    group_by(transcript) %>%
    inner_join(data, by = "transcript") -> DGE_swiss

go %>%
    group_by(transcript) %>%
    inner_join(data, by = "transcript") -> DGE_GO

pfam %>%
    group_by(transcript) %>%
    inner_join(data, by = "transcript") -> DGE_pfam


write.table(DGE_pfam, file=paste0(outpath, "/", file, ".blastx.tsv"), sep="\t", row.names = F, col.names = T)
write.table(DGE_GO, file=paste0(outpath, "/", file, ".go.tsv"), sep="\t", row.names = F, col.names = T)
write.table(DGE_pfam, file=paste0(outpath, "/", file, ".pfam.tsv"), sep="\t", row.names = F, col.names = T)


cat("\nAnnotation tables saved!!!\n")

# The first two columns are: “gene_id” and “transcript_id”, representing predicted genes and their corresponding transcripts, respectively. The columns “Top_BLASTX_hit” and “Top_BLASTP_hit” show the top BlastX and BlastP hit results of homology searches against the NCBI database. BlastX is one of the latest additions to the Trinotate annotation pipeline and compares all six open reading frames (ORF) of the query sequences against the protein database. The RNAMMER column shows information about predicted ribosomal RNA genes discovered in the transcriptome assembly that were predicted by hidden Markov models (HMM). The prot_id, prot_coords and prot_seq columns provide the ID, location and translation of the longest ORFs, respectively. The Pfam column represents the HMMER/PFAM protein domain identification search results. HMMER is used to search databases for homologs of proteins, employing hidden Markov models. The SignalP column shows the presence and location of predicted signal peptides. Similarly, the TmHMM column presents the predicted transmembrane regions. The eggNOG (Evolutionary genealogy of genes: Non-supervised Orthologous Groups) column has the search result of the database of orthologous groups of genes, which are further annotated with functional description lines. Lastly, the gene_ontology column shows the relationship of our data to the Gene Ontology (GO) terms that aim to unify the representation of genes and gene products across all species ----

quit(save = 'no')

DGE_swiss %>%
        select_if(is.numeric) -> Exp

DGE_swiss %>%
        summarise(identity = identity > 50) %>%

        select(transcript, name, uniprot) -> annot

library(RColorBrewer)
library(superheat)

superheat(data[-ncol(data)],
          # retain original order of rows/cols
          pretty.order.rows = TRUE,
          pretty.order.cols = TRUE,
          row.dendrogram = TRUE,
          bottom.label.text.angle = 90,
          row.title = "Differential Expressed",
          column.title = "Samples",
          # change the grid color
          grid.hline.col = "white",
          grid.vline.col = "white",
          left.label.col = "white",
          bottom.label.col = "white",
          heat.pal = brewer.pal(n = 9, name = "Purples"))


cname <- paste0("Sample", 1:ncol(data[-1]))
cname <- paste0("Gene", 1:nrow(data[-1]))

rnames <- aggragate()
aggreg <- aggregate(ftest, list(ftestBestID), sum)

superheat(data[-ncol(data)],
        pretty.order.rows = TRUE,
        label.col = cname)

left.label
bottom.label

genesList <- rownames(data)
genesListgenesList <- go[go$transcript %in% genesList, ]

write.table(go, file=paste0(outpath, "/", "DGE.annot.tsv"), sep="\t", row.names = F, col.names = T)

quit(save= 'no')

gene.df <- clusterProfiler::bitr(genesList, fromType = "GOID",
        toType = c("ENSEMBL", "SYMBOL"),
        OrgDb = org.Hs.eg.db)
head(gene.df)

library(ReactomePA)

x <- enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)




library(org.Hs.eg.db)
library(clusterProfiler)

# go <- as.vector(goProfile$go)

ggo <- groupGO(gene     = goProfile$go,
               OrgDb    = org.Hs.eg.db,
               #ont      = "CC",
               level    = 3,
               readable = TRUE,
               keyType = "GO")




hsGO <- GOSemSim::godata('org.Hs.eg.db', ont="CC", computeIC=FALSE)

go <- as.vector(go$go)
go1 <- sample(go, 20)
go2 <- sample(go, 20)

go <- as.vector(goProfile$go)

table(goProfile$ontology) == max(table(goProfile$ontology))

ontology <- as.data.frame(table(goProfile$ontology) == max(table(goProfile$ontology)))
rownames(ontology) == TRUE

go <- goProfile[goProfile$ontolog == "molecular_function",]
# mgoSim(go1, go2, semData=hsGO, measure="Wang", combine="BMA")

gosim <- GOSemSim::mgoSim(go$go,go$go, semData=hsGO, measure="Wang", combine=NULL)

gosim <- as.data.frame(gosim)

names(gosim) <- go$name

# ========

library(superheat)
superheat(gosim,
          # retain original order of rows/cols
          pretty.order.rows = TRUE,
          pretty.order.cols = TRUE,
          #left.label = "none",
          bottom.label.text.angle = 90,
          row.title = "Sample 1",
          column.title = "Sample 2",
          bottom.label.text.size = 4,
          left.label.text.size = 4
          )

detach(package:org.Hs.eg.db, unload = TRUE)


