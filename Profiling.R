#!/usr/bin/env Rscript

# ========
# Reading Data
# ========

rm(list=ls())

# ========

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
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
  data = read.table(args[2], header=T, com='', row.names=1, check.names=F, sep='\t', stringsAsFactors = FALSE)
  #data$transcript <- rownames(data)
  #data = as.matrix(data)
  
}

# ========================
# Defining variables (some)
# ========================

go <- annot.env$go
blastx <- annot.env$blastx
pfam <- annot.env$pfam
outpath <- getwd()

# quit(save= 'no')

# ===============
# Loading package
# ================
if (!require("dplyr")) { 
        install.packages('dplyr', dep=TRUE, repos='http://cran.us.r-project.org') 
  } else
  if (!require("purrr")) {
      install.packages('purrr', dep=TRUE, repos='http://cran.us.r-project.org')
   } else
  if (!require("tibble")) {
      install.packages('tibble', dep=TRUE, repos='http://cran.us.r-project.org')
   } else
  if (!require("reshape2")) {
      install.packages('reshape2', dep=TRUE, repos='http://cran.us.r-project.org')
   }  else
  if (!require("ggplot2")) {
      install.packages('ggplot2', dep=TRUE, repos='http://cran.us.r-project.org')
   }


# testing ...
#load('Trinotate.xls.RData')
#pathfiles <- "DIFF_EXP/DESeq2.30882.dir/"
#subset <- list.files(path=pathfiles, pattern="UP.subset")
#matrix <- list.files(path=pathfiles, pattern="matrix")
#data = read.table(paste0(pathfiles, matrix[2]), header=T, com='', row.names=1, check.names=F, sep='\t')
#####

data = sample(rownames(data), 10)

blastx %>%
    group_by(transcript) %>%
    inner_join(data, by = "transcript") -> DE_swiss

go %>%
    group_by(transcript) %>%
    inner_join(data, by = "transcript") -> DE_GO

pfam %>%
    group_by(transcript) %>%
    inner_join(data, by = "transcript") -> DE_pfam

head(data)

DE_pfam
DE_GO

quit(save = 'no')



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


