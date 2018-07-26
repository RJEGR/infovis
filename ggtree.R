# powerfull vittage at:
# https://bioconductor.org/packages/devel/bioc/vignettes/ggtree/inst/doc/treeManipulation.html#scale-clade

library("treeio")
library("ggtree")

dir <- c("/Users/cigom/Desktop/viki")
setwd(dir)
# Load tree :::::
#netwick <- args[0]
tree <- treeio::read.newick("bacterial.trim.contigs.ALL.fasta.NJ.tree.2.txt")
tbl <- read.csv("contigs.sina.csv", sep=",", stringsAsFactors = FALSE, header=TRUE)

png(filename = "tree_type.png",  
    width = 20, height = 10, units = "in",
    res = 460)
multiplot(
        ggtree(tree) + geom_treescale(width=0.4),
        ggtree(tree, branch.length = "none") + geom_treescale(width=0.4),
        ggtree(tree, branch.length='none', layout="daylight") + geom_treescale(width=0.4) + geom_nodepoint(),
        ggtree(tree, branch.length='none', layout='circular') + geom_treescale(width=0.4),
        ncol=4,labels = LETTERS[1:4])
dev.off()
# and labels
p <- ggtree(tree) + 
        geom_nodepoint(color="#b5e521", alpha=1/4, size=4) + geom_treescale(width=0.4)


png(filename = "node_type.png",  
    width = 18, height = 10, units = "in",
    res = 460)
multiplot(# add node points
        p + geom_nodepoint(),
        # add tip points
        p + geom_tippoint(),
        # Label the tips
        p + geom_tiplab(),
        ncol=3,labels = LETTERS[1:3])
dev.off()

# other approach
d <- ggtree(tree) + geom_treescale(width=0.4) + geom_nodepoint()

png(filename = "tree.png",  
    width = 16, height = 10, units = "in",
    res = 460)
d %<+% tbl + 
    geom_tiplab(size=4, aes(label=paste0('italic(', label, ')~bolditalic(', Family, ')~', Genus)), parse=TRUE) 
dev.off()

png(filename = "tree.2.png",  
    width = 16, height = 10, units = "in",
    res = 460)
d %<+% tbl + 
    geom_tiplab(size=4, aes(label=paste0('italic(', label, ')~bolditalic(', Family, ')~', Genus)), parse=TRUE) +
    geom_point(aes(color=Family), size=5, alpha=.5) + theme(legend.position="right")

dev.off()



# :: plot aligment
m <- d + geom_tiplab()
# geom_tiplab(align=TRUE, linesize=.5)
fst <- read.fasta("./bacterial.trim.contigs.ALL.fasta.clustalw.2.txt")

png(filename = "tree.3.png",  
    width = 16, height = 10, units = "in",
    res = 460)
msaplot(p=m, fasta=fst, height = 0.4, offset=0.05)
dev.off()
#:::::

png(filename = "tree.nodes.png",  
    width = 20, height = 12, units = "in",
    res = 460)

multiplot(d + geom_treescale(),
        d + geom_text(aes(label=node), hjust=-.3),
        d + geom_tiplab(),
        ncol=3,labels = LETTERS[1:3])
dev.off()

#:::: paint something more
n <- MRCA(tree, tip=c("A7_7", "A7_8"))
n
# [1] 16

png(filename = "tree.further.png",  
    width = 20, height = 12, units = "in",
    res = 460)
multiplot(    
    ggtree(tree) + geom_cladelabel(node = n, label="Random text", color="red"),
    ggtree(tree) + geom_hilight(node = n, "steelblue"),
    ggtree(tree) + geom_hilight(node = n, "steelblue") + geom_cladelabel(node=16, label="Random text", color="red"),
    ncol=3,labels = LETTERS[1:3] )
dev.off()


#:: in python reverse complementary sequence from 
from Bio.Seq import Se
from Bio import SeqIO
for record in SeqIO.parse("../reverse_only.fasta", "fasta"):
    print(record.id)
    print(record.seq.reverse_complement())