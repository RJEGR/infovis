# Solution from https://github.com/benjjneb/dada2/issues/236#issuecomment-422865307
### Qual VS MaxEE Plot
### Scripted edited from Remi Maglione from Kembel Lab


#### Extract qual with FastQC and sed
# fastqc --nogroup (requiered) yourFastqFile_R1.fastq yourFastqFile_R2.fastq

###For R analysis you will need to load the multiplot.R function : 
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

### Function
qualMaxEEplot <- function(fastq.r1, fastq.r2, name1= "Fastq_R1", name2="Fastq_R2") {
  ### Loading dependenies
  .cran_packages <- c('ggplot2', 'fastqcr')
  sapply(c(.cran_packages), require, character.only = TRUE, quietly = TRUE)
  ### Definning function  
  fastqTmp <- function(fileqc) {
    fastq <- data.frame(qc_read(fileqc, modules = "Per base sequence quality", verbose = TRUE)$per_base_sequence_quality)
    fastq.tmp <- rbind(data.frame(R=fastq$Base, 
                                  Q=fastq$Mean, S=c("Mean"), 
                                  E=10^(-fastq$Mean/10), 
                                  A=Reduce('+', 10^(-fastq$Mean/10), accumulate = TRUE)),
                       data.frame(R=fastq$Base, 
                                  Q=fastq$Median, 
                                  S=c("Median"), 
                                  E=10^(-fastq$Median/10), 
                                  A=Reduce('+', 10^(-fastq$Median/10), accumulate = TRUE)),
                       data.frame(R=fastq$Base, 
                                  Q=fastq$Lower.Quartile, 
                                  S=c("Lower.Quartile"), 
                                  E=10^(-fastq$Lower.Quartile/10), 
                                  A=Reduce('+', 10^(-fastq$Lower.Quartile/10), accumulate = TRUE)),
                       data.frame(R=fastq$Base, 
                                  Q=fastq$Upper.Quartile, 
                                  S=c("Upper.Quartile"), 
                                  E=10^(-fastq$Upper.Quartile/10), 
                                  A=Reduce('+', 10^(-fastq$Upper.Quartile/10), accumulate = TRUE)),
                       data.frame(R=fastq$Base, 
                                  Q=fastq$X10th.Percentile, 
                                  S=c("X10th.Percentile"), 
                                  E=10^(-fastq$X10th.Percentile/10), 
                                  A=Reduce('+', 10^(-fastq$X10th.Percentile/10), accumulate = TRUE)),
                       data.frame(R=fastq$Base, 
                                  Q=fastq$X90th.Percentile, 
                                  S=c("X90th.Percentile"), 
                                  E=10^(-fastq$X90th.Percentile/10), 
                                  A=Reduce('+', 10^(-fastq$X90th.Percentile/10), accumulate = TRUE)))
    return(fastq.tmp)
  }
  
  qualPlot <- function(df.tmp) {
    p_r <- ggplot(df.tmp, aes(color=S)) + 
      geom_point(aes(x=R, y=Q), size=1) + 
      labs(x="Reads position", y="Reads Quality")
    return (p_r)
  }
  
  maxEEplot <- function(df.tmp) {
    q_r <- ggplot(df.tmp, aes(color=S)) + 
      geom_point(aes(x=R, y=log10(A)), size=1) + 
      geom_hline(yintercept=log10(2), color = "#de2d26") + 
      geom_hline(yintercept=log10(3), color = "#de2d26") + 
      geom_hline(yintercept=log10(5), color = "#de2d26") + 
      geom_hline(yintercept=log10(7), color = "#de2d26") +
      geom_text(label="MaxEE=2", aes(x=0, y=log10(2), hjust = 0, vjust=0), color="#de2d26") + 
      geom_text(label="MaxEE=3", aes(x=0, y=log10(3), hjust = 0, vjust=0), color="#de2d26") + 
      geom_text(label="MaxEE=5", aes(x=0, y=log10(5), hjust = 0, vjust=0), color="#de2d26") + 
      geom_text(label="MaxEE=7", aes(x=0, y=log10(7), hjust = 0, vjust=0), color="#de2d26") +
      labs(x="Reads position", y="EE = sum(10^(-Q/10)) log10") +
      coord_cartesian(ylim = c(log10(min(df.tmp$A)), log10(max(df.tmp$A)))) + scale_color_brewer(palette = 'Dark2') 
    return (q_r)
  }
  
  ### MAIN  
  p_r1 <- qualPlot(df.tmp = fastqTmp(fastq.r1)) + 
    ggtitle(name1) + scale_color_brewer(palette = 'Dark2') +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  p_r2 <- qualPlot(df.tmp = fastqTmp(fastq.r2)) + ggtitle(name2) + 
    ggtitle(name2) + scale_color_brewer(palette = 'Dark2') +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  q_r1 <- maxEEplot(df.tmp = fastqTmp(fastq.r1))
  q_r2 <- maxEEplot(df.tmp = fastqTmp(fastq.r2))
  
  return(multiplot(p_r1, q_r1, p_r2, q_r2, cols = 2))
}

### Example of usage
# data_path <- '/Users/cigom/metagenomics/COI/run15/'
# fileqc <- list.files(data_path, pattern = 'zip', full.names = TRUE)
# fastq.r1 <- fileqc[1]
# fastq.r2 <- fileqc[2]
# qualMaxEEplot(fastq.r1 = fastq.r1, fastq.r2 = fastq.r2)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
