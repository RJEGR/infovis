# https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0393-0
# Here, we present a novel extension of the SParse InversE Covariance estimation for Ecological ASsociation Inference (SPIEC-EASI) framework that allows statistical inference of cross-domain associations from targeted amplicon sequencing data. By using independent TAS studies of the 16S rRNA gene and Internal Transcribed Spacer (ITS) from the same samples of bacterial and fungal communities, our novel SPIEC-EASI variant allows compositionally robust, simultaneous inference of both within-domain and cross-domain associations

# curated manual here: https://github.com/zdk123/SpiecEasi
# 

library(devtools)
# curl -O https://kingaa.github.io/scripts/mac-fortran.sh
# sh mac-fortran.sh
# install_github("zdk123/SpiecEasi")
library(SpiecEasi)

dataList <- readRDS(paste0(dir, '/phyloseqList.rds'))

# We showed that neighborhood selection (S-E(MB)) outperforms SparCC and CCREPE in terms of recovery of taxon-taxon interactions and global network topology features under almost all tested benchmark scenarios, while covariance selection (S-E(glasso)) performs competitively with and sometimes better than SparCC and CCREPE. https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004226#sec012

# and parallelize

## Parallel multicore ##

pargs2 <- list(rep.num = 20, seed=110620, ncores = 4, thresh = 0.05)

se.res <- spiec.easi(dataList, method='mb', 
                     # lambda.min.ratio=1e-2, 
                     nlambda = 30,
                     sel.criterion = 'stars', 
                     pulsar.select = TRUE, 
                     pulsar.params = pargs2)


tax_table(dataList[[1]])

# dtype <- c(rep(1,ntaxa(hmp216S)), rep(2,ntaxa(hmp2prot)))

g <- adj2igraph(getRefit(se.res))

plot(g, vertex.size=9) # vertex.color=dtype+1, 

library(ggraph)
library(igraph)
library(tidygraph)

# V(g)
# E(g)

# graph <- tbl_graph(nodes = Nodes, edges = Edges, directed = FALSE)

pargs3 <- list(rep.num = 30, seed=110620, ncores = 4, thresh = 0.05)

ps <- prepare_ps(fileNames[1], agg = T)  
  # prune_taxa(taxa_sums(.) > 0, .) %>% 
  # prune_samples(sample_sums(.) > 0, .) %>%
  # transform_sample_counts(., function(x) sqrt(x / sum(x)))

se.mb.ps <- spiec.easi(ps, method='mb', lambda.min.ratio=1e-2,
                          nlambda = 20, sel.criterion = 'stars', 
                       pulsar.select = TRUE, 
                       pulsar.params = pargs3)

# computational cost if not agglomerate

se.gl.ps <- spiec.easi(ps, method='glasso', lambda.min.ratio=1e-2,
                          nlambda=20, pulsar.select = TRUE,
                       pulsar.params = pargs3)

sparcc.ps <- sparcc(otu_table(ps))

hist(abs(sparcc.ps$Cor))

## Define arbitrary threshold for SparCC correlation matrix for the graph

sparcc.graph <- abs(sparcc.ps$Cor) >= 0.1
diag(sparcc.graph) <- 0
library(Matrix)
sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
## Create igraph objects
ig.mb     <- adj2igraph(getRefit(se.mb.ps))
ig.gl     <- adj2igraph(getRefit(se.gl.ps))
ig.sparcc <- adj2igraph(sparcc.graph)

library(igraph)
## set size of vertex proportional to clr-mean
vsize    <- colMeans(clr(otu_table(ps), 1))+6
am.coord <- layout.fruchterman.reingold(ig.mb)

par(mfrow=c(1,3))

plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB")
plot(ig.gl, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc")

