#!/usr/bin/env Rscript

# ========
# Reading RData
# ========
# https://www.r-bloggers.com/safe-loading-of-rdata-files-2/

rm(list=ls())
.LoadToEnvironment <- function(RData, env = new.env()){
  load(RData, env)
  return(env) 
}
annot.env <- .LoadToEnvironment('Trinotate.xls.RData')


# ==============
#  Reading DifExp list
# ==============

args = commandArgs(trailingOnly=TRUE)


if (length(args)==0) {
  stop("!!!\n A matrix of Differential genes Expressed needed, please input in the code as follow: \nRscript Profiling.R diffExpr.*.matrix .\n", call.=FALSE)
} else {
    print('Reading matrix file.')
    data = read.table(arg[1], header=T, com='', row.names=1, check.names=F, sep='\t')
    data = as.matrix(data)
}
#
go <- annot.env$go

# testing ...
load('Trinotate.xls.RData')
pathfiles <- "DIFF_EXP/DESeq2.30882.dir/"
subset <- list.files(path=pathfiles, pattern="UP.subset")
matrix <- list.files(path=pathfiles, pattern="matrix")

data = read.table(paste0(pathfiles, matrix[2]), header=T, com='', row.names=1, check.names=F, sep='\t')
#####

genesList <- rownames(data)
genesListgenesList <- go[go$transcript %in% genesList, ]

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


