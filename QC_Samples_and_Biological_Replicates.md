# QC Samples and Biological Replicates

Once you've performed transcript quantification for each of your biological replicates, it's good to examine the data to ensure that your biological replicates are well correlated, and also to investigate relationships among your samples (Brian Haas, 2018, web:https://github.com/trinityrnaseq/trinityrnaseq/wiki/QC-Samples-and-Biological-Replicates). 

## Batch effects

Batch effects are sub-groups of measure- ments that have qualitatively different behaviour across conditions and are unrelated to the biological or scientific variables in a study. For example, batch effects may occur if a subset of experiments was run on Monday and another set on Tuesday, if two technicians were responsible for different subsets of the experiments or if two different lots of reagents, chips or instruments were used. These effects are not exclusive to high- throughput biology and genomics research, and batch effects also affect low-dimensional molecular measurements, such as northern blots and quantitative PCR. Although batch effects are difficult or impossible to detect in low-dimensional assays, high-throughput technologies provide enough data to detect and even remove them. However, if not properly dealt with, these effects can have a particularly strong and pervasive impact.

One way to quantify the affect of non- biological variables is to **examine the principal components** of the data. Principal components are estimates of the most com- mon patterns that exist across features. For example, if most genes in a microarray study are differentially expressed with respect to cancer status, the first principal compo- nent will be highly correlated with cancer status. Principal components capture both biological and technical variability and, in some cases, principal components can be estimated after the biological variables have been accounted for15. In this case, the prin- cipal components primarily quantify the effects of artefacts on the high-throughput data. Principal components can be com- pared to known variables, such as processing group or time. If the principal components do not correlate with these known variables, there may be an alternative, unmeasured source of batch effects in the data. 

**Tackling the widespread and critical impact of batch effects in high-throughput data (2010) Jeffrey T. Leek, Robert B. Scharpf, Héctor Corrada Bravo, David Simcha, Benjamin Langmead, W. Evan Johnson, Donald Geman, Keith Baggerly and Rafael A. Irizarry. NATURE REVIEWS | GENETICS** 

> If you have replicates that are clear outliers, you might consider removing them from your study as potential confounders. If it's clear that you have a [batch effect](http://www.nature.com/nrg/journal/v11/n10/full/nrg2825.html), you'll want to eliminate the batch effect during your downstream analysis of differential expression.

The trinity utils have the follow code to visualize some intersting patters from the abundance count matrix. Copy and paste the follow code within an script  (ex. `PtR.sh` ) and input both, your abundance count matrix as well as the samples.file used during your analysis. The final code syntax will be `bash PtR.sh count.matrix samples.file`

```bash
#!/bin/bash -ve
PROG=~/Documents/Tools/trinityrnaseq-master/Analysis/DifferentialExpression/PtR
FILE= $1 # input your count matrix
SAMPLE=$2

#if [ -e $FILE.gz ] && [ ! -e $FILE.matrix ]; then
#    gunzip -c $FILE > $FILE.matrix
#fi

# assessing read counts and genes mapped
$PROG --matrix $FILE --samples $SAMPLE --boxplot_log2_dist 1 --output boxplots

# compare replicates
$PROG --matrix $FILE --samples $SAMPLE --compare_replicates --log2  --output rep_compare

# correlation of expression across all samples
$PROG --matrix $FILE --samples $SAMPLE --log2 --sample_cor_matrix --output sample_cor

# PCA analysis
$PROG --matrix $FILE --samples $SAMPLE --log2 --prin_comp 4 --output prin_comp

# Most variably expressed genes
$PROG --matrix $FILE --samples $SAMPLE --log2 --top_variable_genes 1000 --var_gene_method anova --heatmap --center_rows --output anovaTop1k

# combine analyses
$PROG --matrix $FILE --samples $SAMPLE --log2 \
      --top_variable_genes 1000 --var_gene_method anova --output anovaTop1k --heatmap --prin_comp 3  --add_prin_comp_heatmaps 30 \
      --center_rows --output top1kvarPC3Gene30
```

> PCA is a method for compressing a lot of data into something that capture the essence of the original data. 
>
>  The goal of PCA is to draw a graph that shows how the samples are related (or not related) to each other
>
> Técnicamente, el ACP busca la proyección según la cual los datos queden mejor representados en términos de mínimos cuadrados. Esta convierte un conjunto de observaciones de variables posiblemente correlacionadas en un conjunto de valores de variables sin correlación lineal llamadas **componentes principales** (wikipedia)

## An introduction to dimension (By Josh Starmer)

RNA-seq results often contain a PCA (Principal Component Analysis) or MDS plot. Usually we use these graphs to verify that the control samples cluster together. However, there’s a lot more going on, and if you are willing to dive in, you can extract a lot more information from these plots. The good news is that PCA only sounds complicated. Conceptually, it’s actually quite simple (StatQuest, 2018). 

* PCA can take 4 or more genes (or measurements) ie. four or more dimensions of data and make a 2-D PCA plot. 
* PCA can tell us which gene (or variable) is the most valuable for clustering the data.
* PCA can tell us how accurate the 2-D graph is.
* PCA is a linear combination variables

Principal Component Analysis, is one of the most useful data analysis and machine learning methods out there. It can be used to identify patterns in highly complex datasets and it can tell you what variables in your data are the most important. Lastly, it can tell you how accurate your new understanding of the data actually is.

Video: [StatQuest: Principal Component Analysis (PCA) clearly explained (2015)](https://www.youtube.com/watch?v=_UVHneBUBW0&t=621s)

### Principal component example

In this chunk code from josh Starmer (StatQuest)  we will create a matrix with colnames (the samples name of wt and ko factor) and rownames (the isoform / gene reported transcript ) 

```R
## columns are individual samples (i.e. cells)
## rows are measurements taken for all the samples (i.e. genes)

data.matrix <- matrix(nrow=100, ncol=10)


colnames(data.matrix) <- c(
  paste("wt", 1:5, sep=""),
  paste("ko", 1:5, sep=""))


rownames(data.matrix) <- paste("isoform_s", 1:100, sep="")
```

Then let's reproduce a random possion distribution of values than been used as counts for each of the isoforms/genes

```R


for (i in 1:100) {
  wt.values <- rpois(5, lambda=sample(x=10:1000, size=1))
  ko.values <- rpois(5, lambda=sample(x=10:1000, size=1))
  
  data.matrix[i,] <- c(wt.values, ko.values)
}

head(data.matrix)


```

The head of the data.matrix looks like:

|            |  wt1 |  wt2 |  wt3 |  wt4 |  wt5 |  ko1 |  ko2 |  ko3 |  ko4 |  ko5 |
| :--------- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| isoform_s1 |  152 |  145 |  188 |  176 |  154 |  381 |  407 |  383 |  375 |  402 |
| isoform_s2 |  383 |  416 |  379 |  409 |  403 |  769 |  681 |  769 |  782 |  740 |
| isoform_s3 |  848 |  842 |  843 |  920 |  876 |  184 |  191 |  180 |  216 |  217 |
| isoform_s4 |  798 |  718 |  777 |  744 |  738 |  123 |  124 |  118 |  112 |   73 |
| isoform_s5 |  565 |  512 |  534 |  504 |  513 |  841 |  836 |  835 |  805 |  776 |
| isoform_s6 |  182 |  170 |  176 |  180 |  160 |  941 |  874 |  916 |  875 |  866 |

The final goal from this chuck is to dra a graph that shows how the samples are related (or not related) to each other. Using the `prcomp()` function in r. This function retur three things: 

1. x (where x contain the PCs for drawing a graph. Since there are 10 samples, there are 10 PCs).
2. Standar Deviations (sdev).
3. number of rotation.

By default prcomp expects the samples to be rows and genes to be columns. Let's transpose the data.matrix yet and plot the Principal component using the grammar of graphics with ggplot.

```R

pca <- prcomp(t(data.matrix), scale=TRUE)

library(ggplot2)

pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("My PCA Graph")



```

> The first PC accounts for the most variation in the original data (the gene expression across all the 10 samples). The 2nd PC accounts for the second most variation and so on.
>
> To plot a 2-D PCs graph, we usually use the fist PCs. sometimes we use the PC1 and PC3

To get sense of how meaninfull these clusters are, let's see how much variation in the original data PC1 accounts for. We use the square of sdev, which stand for _standar deviation_ to calculate how much variation in the orignal data each principal component accounts for. Since the percent of variation that each PC accounts for, is way more intersting that the actual value, we calculate the percentages. Finally plot the results.

```r

# 1
pca.var <- pca$sdev^2

# 2
pca.var.per <- round(pca.var/sum(pca.var)*100,1 ) 

# 3. and plot
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")


```

Also we can figure out wich isoforms (variables) contribute to the variation of the PCs. 

```R
## get the name of the sample (cell) with the highest pc1 value
rownames(pca$x)[order(pca$x[ ,1], decreasing=TRUE)[1]]

## get the name of the top 10 measurements (genes) that contribute
## most to pc1.
loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_10_genes <- names(gene_score_ranked[1:10])

top_10_genes ## show the names of the top 10 genes




```



## IGV

Load in the 'Trinity.fasta' file as a 'genome' via the IGV 'Genomes'->'Load Genome from File' menu.

Load in the 'bowtie2.coordSorted.bam' file via the IGV 'File'->'Load from File' menu. (Again, be sure this is the coordinate-sorted version and NOT the name-sorted version.)

> From your earlier expression data (counts matrix), try finding a highly expressed transcript and view alignments to it within IGV. 
>
> ```bash
> grep
> ```



