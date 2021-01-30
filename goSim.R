# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GOSim")

library(GOSim)

getTermSim(c("GO:0007166","GO:0007267","GO:0007584","GO:0007165","GO:0007186"), method="Resnik",verbose=FALSE)


G = getGOGraph(c("GO:0007166","GO:0007267"))

library(igraph)

G2 = igraph.from.graphNEL(G)



plot(G2, vertex.label=V(G2)$name)

getTermSim(c("GO:0007166","GO:0007267"),method="CoutoResnik",verbose=FALSE)
