#!/usr/bin/env Rscript
# Rscript --vanilla annotation.R Trinotate.xls;

# Loadding package:

if (!require('devtools')) {
    install.packages("devtools", dep=TRUE, repos='http://cran.us.r-project.org')
} else   
    if (!require('trinotateR')) {
        devtools::install_github("cstubben/trinotateR")
    } else
    if (!require("DT")) { 
        install.packages('DT', dep=TRUE, repos='http://cran.us.r-project.org') 
    } else
    if (!require("htmlwidgets")) {
        install.packages('htmlwidgets', dep=TRUE, repos='http://cran.us.r-project.org') 
    } else
    if (!require("ggpubr")) {
        install.packages('ggpubr', dep=TRUE, repos='http://cran.us.r-project.org')
     } else
    if (!require("ggplot2")) {
        install.packages('ggplot2', dep=TRUE, repos='http://cran.us.r-project.org')
    }
    

args = commandArgs()
file = args[7]

paste0("Using ", file, " as input file")

# ==== Create a modified function
split_blast <- function (x, hit = "sprot_Top_BLASTX_hit")
{
    y <- x[!is.na(get(hit)), .(get(hit), gene_id, transcript_id,
        prot_id)]
    z <- strsplit(y$V1, "`")
    n <- sapply(z, length)
    z <- strsplit(unlist(z), "\\^")
    if (any(sapply(z, "[", 1) != sapply(z, "[", 2)))
        print("WARNING: check different values in columns 1 and 2")
    NAME <- gsub("^RecName: Full=", "", sapply(z, "[", 6))
    NAME <- gsub("SubName: Full=", "", NAME)
    NAME <- gsub(";$", "", NAME)
    NAME <- gsub(" \\{[^}]+}", "", NAME)
    x1 <- data.frame(gene = rep(y$gene_id, n), transcript = rep(y$transcript_id,
        n), protein = rep(gsub(".*\\|", "", y$prot_id), n), uniprot = sapply(z,
        "[", 1), align = sapply(z, "[", 3), identity = as.numeric(gsub("%ID",
        "", sapply(z, "[", 4))), evalue = as.numeric(gsub("E:",
        "", sapply(z, "[", 5))), name = NAME, lineage = sapply(z,
        "[", 7), domain = gsub("; .*", "", sapply(z, "[", 7)),
        genus = gsub(".*; ", "", sapply(z, "[", 7)), stringsAsFactors = FALSE)
    message(nrow(x1), " ", hit, " annotations")
    data.table(x1)
}
# ====

library(trinotateR)

x <- read_trinotate(file)
summary_trinotate(x)

pfam <- split_pfam(x)
spfam <- summary_pfam(pfam)

go <- split_GO(x)
gos <- summary_GO(go)

blastx <- split_blast(x, "sprot_Top_BLASTX_hit")
sblastx <- summary_blast(blastx)

blastp <- split_blast(x, "sprot_Top_BLASTP_hit")
sblastp <- summary_blast(blastp)

library(DT)

z <- data.frame(spfam)
z$pfam <- paste0('<a href="http://pfam.xfam.org/family/', z$pfam, '">', z$pfam,  '</a>')

widget <- datatable(
  z, 
  escape=1,
  extensions = 'Buttons', options = list(
    pageLength = 25,  
    dom = 'Bfrtip',
    buttons = 
      list('copy', 'print', list(
        extend = 'collection',
        buttons = c('csv', 'excel'),
        text = 'Download'
      ))
    
  )
)

htmlwidgets::saveWidget(widget, paste0(file, ".pfam.html"))


#
data(cogs)
download.file("http://eggnogdb.embl.de/download/eggnog_4.5/data/NOG/NOG.annotations.tsv.gz", "NOG.annotations.tsv.gz")
system("gunzip NOG.annotations.tsv.gz")
egg <- read.table("NOG.annotations.tsv", sep="\t", stringsAsFactors=FALSE, quote="")
names(egg) <- c("db", "nog", "proteins", "species", "class", "description")



plot_NOGs(x, "transcript_id")

system("rm NOG.annotations.tsv")
# or =======
#y <- unique(x[!is.na(eggnog), .(gene_id, eggnog)])
#xlabel <- "Number of genes"

#nogs <- gsub("(.*)\\^.*", "\\1", y$eggnog)
#n <- match(nogs, eggnog$nog)
#y <- table(unlist(strsplit(eggnog$class[n], "")))
#y1 <- rev(y[order(match(names(y), cogs$code))])
#y1 <- y1[-1]
#n <- match(names(y1), cogs$code)

#op <- par(mar = c(4, 25, 1, 1), mgp = c(2, 0.5, 0))
#barplot(y1, col = cogs$clrs[n], horiz = TRUE, names.arg = paste(cogs$code[n],
#        cogs$name[n], sep = ". "), las = 1, xlab = xlabel)
#par(op)

# ====
library(ggpubr)

plotgos <- head(gos[order(-gos$transcripts),], 80)

bar <- ggbarplot(plotgos, "name", "transcripts",
          fill = "ontology", 
          color = "ontology",
          xlab = "Ontology",
          ylab = "Number of transcripts",
   palette = c("#00AFBB", "#E7B800", "#FC4E07"))

#png(filename = paste0(file, 'plot_NOGs.png'), height = 28,width = 36,res = 300,units = "cm")
bar + theme(axis.text.x = element_text(angle = 75, hjust = 1, size = 7))
#dev.off()

# ====== blastp or blastp
genus <- data.table(table(sblastp$genus))
names(genus) <- c("Genus", "Number")
genus <- genus[order(-genus$Number), ]
genus[(Number <= 10 ), Genus := "Others"]

library(ggplot2)

ggplot(genus, aes(x = reorder(Genus, -Number), y = Number)) +
 geom_bar(stat="identity", fill = "steelblue") +
 coord_flip() + xlab("Genus") + ylab("Number of ORF") +
 theme_minimal()

# ========
genus <- NULL
genus <- data.table(table(sblastx$genus))
names(genus) <- c("Genus", "Number")
genus <- genus[order(-genus$Number), ]
genus[(Number <= 10 ), Genus := "Others"]

library(ggplot2)

ggplot(genus, aes(x = reorder(Genus, -Number), y = Number)) +
 geom_bar(stat="identity", fill = "steelblue") +
 coord_flip() + xlab("Genus") + ylab("Number of transcripts") +
 theme_minimal()
 #geom_text(data=subset(genus, Genus=="Others"), aes(label=sum(Number)), hjust = 0, vjust = 1)

paste0("Summary done from ", file, " file")

quit(save = "no")