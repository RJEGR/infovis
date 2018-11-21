#!/usr/bin/env Rscript

# Use:  Rscript --vanilla annotation.R Trinotate.xls;

# ================
# Defining outputs
# ================
outpath <- getwd()

# ===============
# Loading package:
# ===============
lib.loc = "/home/rgomez/R/x86_64-pc-linux-gnu-library/3.5"

.cran_packages <- c("DT", "htmlwidgets", "ggpubr", "ggplot2", "dplyr")
.bioc_packages <- c("devtools", "trinotateR")

.inst <- .cran_packages %in% installed.packages(lib.loc = lib.loc)
if(any(!.inst)) {
   install.packages(.cran_packages[!.inst], dep=TRUE, repos='http://cran.us.r-project.org')
}
.inst <- .bioc_packages %in% installed.packages(lib.loc = lib.loc)

#if(any(!.inst)) {
#    if (!requireNamespace("BiocManager", quietly = TRUE))
#        install.packages("BiocManager")
#    BiocManager::install(.bioc_packages[!.inst], ask = F)
#}

if (any(!.inst)) {   
    if (!require('devtools', lib.loc = lib.loc)) {
        install.packages("devtools", dep=TRUE, repos='http://cran.us.r-project.org')
    } else 
    devtools::install_github("cstubben/trinotateR")
    require('trinotateR')
    }

# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, lib.loc = lib.loc, character.only = TRUE)


# Loadding data
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  file <- "Trinotate.xls"
} else 
  file = args[1]

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

#
data(cogs)
download.file("http://eggnogdb.embl.de/download/eggnog_4.5/data/NOG/NOG.annotations.tsv.gz", "NOG.annotations.tsv.gz")
system("gunzip NOG.annotations.tsv.gz")
egg <- read.table("NOG.annotations.tsv", sep="\t", stringsAsFactors=FALSE, quote="")
names(egg) <- c("db", "nog", "proteins", "species", "class", "description")



plot_NOGs(x, "transcript_id")

system("rm NOG.annotations.tsv")

# ====
library(ggpubr)

plotgos <- head(gos[order(-gos$transcripts),], 80)

bar <- ggbarplot(plotgos, "name", "transcripts",
          fill = "ontology", 
          color = "ontology",
          xlab = "Ontology",
          ylab = "Number of transcripts",
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in dscending order
          sort.by.groups = TRUE)
   #palette = c("#00AFBB", "#E7B800", "#FC4E07"))


#png(filename = paste0(file, 'plot_NOGs.png'), height = 28,width = 36,res = 300,units = "cm")
bar + theme(axis.text.x = element_text(angle = 75, hjust = 1, size = 7))
#dev.off()

# ====== blastp or blastp
genus <- data.table(table(sblastp$genus))
names(genus) <- c("Genus", "Number")
genus <- genus[order(-genus$Number), ]
genus[(Number <= 10 ), Genus := "Others"]

#library(ggplot2)
#ggplot(genus, aes(x = reorder(Genus, -Number), y = Number)) +
# geom_bar(stat="identity", fill = "steelblue") +
# coord_flip() + xlab("Genus") + ylab("Number of ORF") +
# theme_minimal()

# ========
xgenus <- data.table(table(sblastx$genus))
names(xgenus) <- c("Genus", "Number")
xgenus <- xgenus[order(-xgenus$Number), ]
xgenus[(Number <= 10 ), Genus := "Others"]


xgenus$Type <-"Transcript"
genus$Type <-"ORF"

plot <- rbind(genus, xgenus)

ggbarplot(plot, x = "Genus", y = "Number",
          fill = "Type",               # change fill color by cyl
          color = "Type",            # Set bar border colors to white
          palette = "Paired",            # jco journal color palett. see ?ggpar
          group = "Type",
          sort.by.groups = TRUE,      # Sort inside each group
          x.text.angle = 90,           # Rotate vertically x axis texts
          ylab = "Number of annotations",
          rotate = TRUE,
          ggtheme = theme_minimal()
          #facet.by = "Type"
          )
 
 #geom_text(data=subset(genus, Genus=="Others"), aes(label=sum(Number)), hjust = 0, vjust = 1)


x <- blastx$identity
y <- blastp$identity


h1 <- hist(x, breaks=100, col=rgb(1,0,0,1/4), 
                xlab="Identity of aligment (Swissprot - blastp hits)", 
                main="Distribution of identity")

# Add a Normal Curve (Thanks to Peter Dalgaard)
xfit<-seq(min(x),max(x),length=100) 
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
yfit <- yfit*diff(h1$mids[1:2])*length(x) 

# plot( h1, col=rgb(1,0,0,1/4), xlab="Identity of aligment (Swissprot - blastp hits)", main="Distribution of identity")
# plot( h2, col=rgb(1,0,0,1/6), add=T) 
lines(xfit, yfit, col="blue", lwd=2)



z <- data.frame(spfam)
z$pfam <- paste0('<a href="http://pfam.xfam.org/family/', z$pfam, '">', z$pfam,  '</a>')

# ================
# Save image first
save.image(file = paste0(file, ".RData"))

# ======================
# Save annotation outputs
# ======================

write.table(blastx, file=paste0(outpath, "/", file, ".blastx.tsv"), sep="\t", row.names = F, col.names = T)
write.table(go, file=paste0(outpath, "/", file, ".go.tsv"), sep="\t", row.names = F, col.names = T)


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


paste0("\nSummary done from ", file, " file")

quit(save = "no")




# =============
# ==== subset list of differential gene list... next release
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop()
} else 
  data = args[2]
# checar los scritps de diffExp de trinity para saber como hacer el brake stop

read.csv(data, sep = "\t", header = FALSE) %>%
            as.tibble() %>%
            rename(transcript = X1) -> DiffEx

# DiffEx <- data.frame(transcript = sample(blastx$transcript, 40))



blastx %>%
    group_by(transcript) %>%
    inner_join(DiffEx, by = "transcript") -> DE_swiss

go %>%
    group_by(transcript) %>%
    inner_join(DiffEx, by = "transcript") -> DE_GO

pfam %>%
    group_by(transcript) %>%
    inner_join(DiffEx, by = "transcript") -> DE_pfam

#
#

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
    } else
    if (!require("dplyr")) {
        install.packages('dplyr', dep=TRUE, repos='http://cran.us.r-project.org')
    }