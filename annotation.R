#!/usr/bin/env Rscript

# Use:  Rscript --vanilla annotation.R Trinotate.xls;

# ================
# Defining outputs
# ================
outpath <- getwd()

# ===============
# Loading package:
# ===============

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
    
# Loadding data
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  file <- "Trinotate.xls"
} else 
  file = args[1]

cat("\nUsing", file, "as input file\n")

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

# ==== Read trinotate format

x <- read_trinotate(file)
summary_trinotate(x)

# Summary data usage

summary.tbl <- summary_trinotate(x)
summary.tbl <- summary.tbl[summary.tbl[1] > 0,][1]
ggdata <- data.frame(x=rownames(summary.tbl), y=summary.tbl$unique)

ggdata %>%
 arrange(desc(y)) %>%
 mutate(x=factor(x,x)) -> ggdata

# step 2 

data.usage <- c("transcript_id","sprot_Top_BLASTX_hit", "sprot_Top_BLASTP_hit", "Pfam", "gene_ontology_blast", "eggnog")

p <- ggplot(ggdata, aes(x=x, y=y)) +
  geom_segment( aes(x=x, xend=x, y=0, yend=y ), 
                color = ifelse(ggdata$x %in% data.usage, "orange", "grey"), 
                size = ifelse(ggdata$x %in% data.usage, 1.3, 0.7) ) +
  geom_point( color = ifelse(ggdata$x %in% data.usage, "orange", "grey"), 
              size=ifelse(ggdata$x %in% data.usage, 5, 2) ) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="none",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  xlab("") +
  ylab("Number of Sequences") +
  ggtitle("Unigenes Data distribution \nData usage through this analysis is colored with orange")


print(p)

# Split data usage

pfam <- split_pfam(x)
spfam <- summary_pfam(pfam)

go <- split_GO(x)
gos <- summary_GO(go)

blastx <- split_blast(x, "sprot_Top_BLASTX_hit")
sblastx <- summary_blast(blastx)

blastp <- split_blast(x, "sprot_Top_BLASTP_hit")
sblastp <- summary_blast(blastp)


cat("\nRibosomal unit reported in the assembly:\n")
na.exclude(top_table(x, "RNAMMER", n = 10))

#
data(cogs)
download.file("https://raw.githubusercontent.com/RJEGR/infovis/master/NOG.annotations.tsv", "NOG.annotations.tsv")
egg <- read.table("NOG.annotations.tsv", sep="\t", stringsAsFactors=FALSE, quote="")
names(egg) <- c("db", "nog", "proteins", "species", "class", "description")



plot_NOGs(x, "transcript_id")

system("rm NOG.annotations.tsv")

# ====
library(ggpubr)
plotgos <- head(gos[order(-gos$transcripts),], 100)

bar <- ggbarplot(plotgos, "name", "transcripts",
          fill = "ontology", 
          color = "ontology",
          xlab = "Name",
          ylab = "Number of transcripts",
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in dscending order
          sort.by.groups = TRUE,
          title = "Highest ontologies (top 100)")
          
   #palette = c("#00AFBB", "#E7B800", "#FC4E07"))


#png(filename = paste0(file, 'plot_NOGs.png'), height = 28,width = 36,res = 300,units = "cm")
bar + theme(axis.text.x = element_text(angle = 75, hjust = 1, size = 10))
#dev.off()

# ====== blastp or blastp
genus <- data.table(table(sblastp$genus))
names(genus) <- c("Genus", "Number")
genus <- genus[order(-genus$Number), ]
genus[(Number <= round(mean(genus$Number))), Genus := "Others"]

#length(table(genus$Genus))
#library(ggplot2)
#ggplot(genus, aes(x = reorder(Genus, -Number), y = Number)) +
# geom_bar(stat="identity", fill = "steelblue") +
# coord_flip() + xlab("Genus") + ylab("Number of ORF") +
# theme_minimal()

# ========
xgenus <- data.table(table(sblastx$genus))
names(xgenus) <- c("Genus", "Number")
xgenus <- xgenus[order(-xgenus$Number), ]
xgenus[(Number <= round(mean(xgenus$Number))), Genus := "Others"]
#length(table(xgenus$Genus))

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
          ) + theme(axis.text.y = element_text(hjust = 1, size = 7))
 
 #geom_text(data=subset(genus, Genus=="Others"), aes(label=sum(Number)), hjust = 0, vjust = 1)


xy <- data.frame(Identity=blastx$identity, Type="Transcript")
y <- data.frame(Identity=blastp$identity, Type="ORF")

data <- rbind(xy,y)

ggplot(data, aes(Identity, fill = Type)) + 
        geom_histogram(bins = 100, alpha = 0.7, aes(y = ..density..), position = 'identity') + 
        theme_classic() +
        scale_fill_brewer(direction = -1, palette = "Paired") + 
        stat_function(fun=dnorm,
                     color="red",
                     args=list(mean=mean(data$Identity, na.rm = TRUE), 
                              sd=sd(data$Identity, na.rm = TRUE))) +  
        scale_x_continuous("Identity of the aligment")
#geom_histogram(binwidth = 2.5)+
#      geom_density(aes(y=2.5 * ..count..)) <--- para emparejar ambos, histograma y densidad
# x <- blastx$identity
# y <- blastp$identity

# h1 <- hist(x, breaks=100, col=rgb(1,0,0,1/4), 
#                xlab="Identity of aligment (Swissprot - blastp hits)", 
#                main="Distribution of identity")

# Add a Normal Curve (Thanks to Peter Dalgaard)
# xfit<-seq(min(x),max(x),length=100) 
# yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
# yfit <- yfit*diff(h1$mids[1:2])*length(x) 

# plot( h1, col=rgb(1,0,0,1/4), xlab="Identity of aligment (Swissprot - blastp hits)", main="Distribution of identity")
# plot( h2, col=rgb(1,0,0,1/6), add=T) 
#lines(xfit, yfit, col="blue", lwd=2)

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


cat("\nSummary done from ", file, " file\n")

quit(save = "no")