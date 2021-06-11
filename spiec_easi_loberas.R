# https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0393-0
# Here, we present a novel extension of the SParse InversE Covariance estimation for Ecological ASsociation Inference (SPIEC-EASI) framework that allows statistical inference of cross-domain associations from targeted amplicon sequencing data. By using independent TAS studies of the 16S rRNA gene and Internal Transcribed Spacer (ITS) from the same samples of bacterial and fungal communities, our novel SPIEC-EASI variant allows compositionally robust, simultaneous inference of both within-domain and cross-domain associations

# curated manual here: https://github.com/zdk123/SpiecEasi

library(devtools)
# curl -O https://kingaa.github.io/scripts/mac-fortran.sh
# sh mac-fortran.sh
# install_github("zdk123/SpiecEasi")
library(SpiecEasi)

dataList <- readRDS(paste0(dir, '/phyloseqList.rds'))

# We showed that neighborhood selection (S-E(MB)) outperforms SparCC and CCREPE in terms of recovery of taxon-taxon interactions and global network topology features under almost all tested benchmark scenarios, while covariance selection (S-E(glasso)) performs competitively with and sometimes better than SparCC and CCREPE. https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004226#sec012

# and parallelize

## Parallel multicore ##
pargs2 <- list(rep.num = 30, seed=110620, ncores = 4, thresh = 0.05)

time <- system.time(
  se.res <- spiec.easi(dataList, method='mb', 
                       lambda.min.ratio=1e-3, nlambda=30,
                    sel.criterion = 'bstars', 
                    pulsar.select = TRUE, 
                    pulsar.params = pargs2)
)




dtype <- c(rep(1,ntaxa(hmp216S)), rep(2,ntaxa(hmp2prot)))

graphdf <- adj2igraph(getRefit(se.res))

plot(graphdf, vertex.color=dtype+1, vertex.size=9)

# data(hmp2)
# se.hmp2 <- spiec.easi(list(hmp216S, hmp2prot), method='mb', nlambda=40,
#                       lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05))
# 
# dtype <- c(rep(1,ntaxa(hmp216S)), rep(2,ntaxa(hmp2prot)))
# 
# graphdf <- adj2igraph(getRefit(se.hmp2))
# 
# plot(graphdf, vertex.color=dtype+1, vertex.size=9)
# 
# c(taxa_names(hmp216S), taxa_names(hmp2prot))
# 
# data('amgut2.filt.phy')
# se.mb.amgut2 <- spiec.easi(amgut2.filt.phy, method='mb', lambda.min.ratio=1e-2,
#                            nlambda=20, pulsar.params=list(rep.num=50))
# ig2.mb <- adj2igraph(getRefit(se.mb.amgut2),  vertex.attr=list(name=taxa_names(amgut2.filt.phy)))
# plot_network(ig2.mb, amgut2.filt.phy, type='taxa', color="Rank3")
