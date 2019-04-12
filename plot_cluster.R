#!/usr/bin/env Rscript
# Rscript --vanilla plot_cluster.R Cluster_Anotación-1a.txt Cluster_Anotación-1a.png
# Mayo 2019
# Ricardo Gomez-Reyes

rm(list = ls())

# ==============
## Checking and Load packages ----
# ==============
.cran_packages <- c('ggplot2', "reshape2")

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst], dep=TRUE, repos='http://cran.us.r-project.org')
}

# Load packages into session, and print package version
sapply(c(.cran_packages), require, character.only = TRUE)

args = commandArgs(trailingOnly=TRUE)
file_in = args[1]
file_out = args[2]

path = getwd()
out_path = paste0(path, "/", "clusters_fig")
system(command = paste0("mkdir -p ", out_path), intern = F)

x <- read.table(file_in, header = TRUE)

levels <- c(names(x[-1]))

xplot <- melt(t(x[-1]))

xplot$Var1 <- factor(xplot$Var1, levels = levels)

pcp_cl <- ggplot(xplot,   aes(Var1, value, group = Var2) )         
jit<- position_jitter(width = .08, height = .08)

#Ok first graph the cluster means.

pcp_cl + stat_summary(fun.y = mean,   geom = "line")

#Then we produce a colourful but uninformative parallel coordinates 
#plot with a bit of alpha blending and jitter.
scale_x_dis <- c("WB-24°C",	"WB-30°C",	"OL-24°C",	"OL-30°C",	"OG-24°C",	"OG-30°C")

save_plot <- pcp_cl + geom_line(position = jit, alpha = 1/5) + 
             stat_summary(fun.y = mean,   geom = "line", colour = "red", aes(group = 1)) +
             theme_classic()  +
             labs(x = "Treatment", y= "Log2FC",
                  title = paste0(nrow(x), " Transcripts")) +
              scale_x_discrete(labels=scale_x_dis)
# transcript_id	c("WB-24°C",	"WB-30°C",	"OL-24°C",	"OL-30°C",	"OG-24°C",	"OG-30°C")

ggsave(paste0(file_out, ".png"), plot=save_plot, path=out_path, width = 10, height = 7, units = "in")

cat('/n Cluster plot Done! /n')
# https://adventuresinr.wordpress.com/2010/10/22/displaying-k-means-results/