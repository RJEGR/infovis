# Genorm, normfinder y bestkeep

# .cran_packages <- c("wTO", "CoDiNA") # "tidyverse"
.bioc_packages <- c("HTqPCR")

# .inst <- .cran_packages %in% installed.packages()
# if(any(!.inst)) {
#   install.packages(.cran_packages[!.inst], dep=TRUE, repos='http://cran.us.r-project.org')
# }
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(.bioc_packages[!.inst], ask = F)
}

require(HTqPCR)

data(qPCRraw)

slotNames(qPCRraw)
pData(qPCRraw)
phenoData(qPCRraw)
featureData(qPCRraw)
head(fData(qPCRraw))

g <- featureNames(qPCRraw)[1:10]

plotCtReps(qPCRraw, card = 2, percent = 20)
plotCtCard(qPCRraw, plot = "class", well.size = 2.6)
plotCtCard(qPCRraw, col.range = c(10, 35), well.size = 2.6)
plotCtOverview(qPCRraw, genes = g, xlim = c(0, 50), conf.int = TRUE, ylim = c(0, 55))

