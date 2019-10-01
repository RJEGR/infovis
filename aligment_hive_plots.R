library("igraph")
library("plyr")
library("HiveR")
library("RColorBrewer")

rm(list = ls())
url_file <- c("http://www.vesnam.com/Rblog/wp-content/uploads/2013/07/lesmis.txt")


# Read a data set. 
# Data format: dataframe with 3 variables; variables 1 & 2 correspond to interactions; variable 3 corresponds to the weight of interaction

dataSet <- read.table(url_file, header = FALSE, sep = "\t")


# Create a graph. Use simplify to ensure that there are no duplicated edges or self loops
gD <- simplify(graph.data.frame(dataSet, directed=FALSE))

# Print number of nodes and edges
vcount(gD)
ecount(gD)

degAll <- degree(gD, v = V(gD), mode = "all")


# Calculate betweenness for all nodes
betAll <- betweenness(gD, v = V(gD), directed = FALSE) / (((vcount(gD) - 1) * (vcount(gD)-2)) / 2)
betAll.norm <- (betAll - min(betAll))/(max(betAll) - min(betAll))
rm(betAll)

# Calculate Dice similarities between all pairs of nodes
dsAll <- similarity.dice(gD, vids = V(gD), mode = "all")

############################################################################################
# Add new node/edge attributes based on the calculated node properties/similarities

gD <- set.vertex.attribute(gD, "degree", index = V(gD), value = degAll)
gD <- set.vertex.attribute(gD, "betweenness", index = V(gD), value = betAll.norm)

# Check the attributes
# summary(gD)

F1 <- function(x) {data.frame(V4 = dsAll[which(V(gD)$name == as.character(x$V1)), which(V(gD)$name == as.character(x$V2))])}
dataSet.ext <- ddply(dataSet, .variables=c("V1", "V2", "V3"), function(x) data.frame(F1(x)))

gD <- set.edge.attribute(gD, "weight", index = E(gD), value = 0)
gD <- set.edge.attribute(gD, "similarity", index = E(gD), value = 0)

# The order of interactions in gD is not the same as it is in dataSet or as it is in the edge list,
# and for that reason these values cannot be assigned directly

E(gD)[as.character(dataSet.ext$V1) %--% as.character(dataSet.ext$V2)]$weight <- as.numeric(dataSet.ext$V3)
E(gD)[as.character(dataSet.ext$V1) %--% as.character(dataSet.ext$V2)]$similarity <- as.numeric(dataSet.ext$V4)

# 
#install.packages('rgexf')
library("rgexf")

# Create a dataframe nodes: 1st column - node ID, 2nd column -node name
nodes_df <- data.frame(ID = c(1:vcount(gD)), NAME = V(gD)$name)
# Create a dataframe edges: 1st column - source node ID, 2nd column -target node ID
edges_df <- as.data.frame(get.edges(gD, c(1:ecount(gD))))

