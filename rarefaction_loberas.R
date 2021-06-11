
ranks <- c('Kingdom',  'Phylum',  'Class',  'Order', 'Family', 'Genus', 'Species', 'pplacer_sp')
sampleName <- c('Cantiles', 'Coloradito', 'Granito', 'Machos', 'Partido', 'Rasito')
dir <- '~/metagenomics/Loberas_MG/'

fileNames <- list.files(dir, pattern = "xls", full.names = T)

    
prepare_ps <- function(filename) {
  
  ranks <- c('Kingdom',  'Phylum',  'Class',  'Order', 'Family', 'Genus', 'Species')
  
  readXLS <- function(x) {
    group <- paste(sapply(strsplit(basename(x), "_"), `[`, 1),collapse = '_')
    readxl::read_xlsx(x, na = 'NA') %>% mutate(group = group)
  }

  features <- readXLS(filename)
  
  names(features)[9] <- 'pplacer_sp'
  
  features %>% 
    mutate(Relationship = recode_factor(Relationship, P = 'Pathogenic', 
                                        NC = 'Inconsistent', C = 'non-pathogenic')) %>%
    mutate_all(., funs(str_replace_all(., c("Bacteroides tectu$" = "Bacteroides tectus")))) %>% 
    mutate_at(sampleName, as.double) -> features
  
  dat <- features %>% select_at(sampleName) %>% data.frame(row.names = 1:nrow(features))
  
  features %>% select(Relationship, pplacer_sp) -> pplacer
  
   # hay problemas con la taxonomia que francesco curo, por tanto tener cuidado al usar datos aglomerados por taxonomia
  
  features %>% 
    select_at(all_of(ranks)) %>%
    mutate(id = 1:nrow(features)) %>%
    pivot_longer(cols = ranks) %>% fill(value) %>% 
    pivot_wider(names_from = name) %>%
    cbind(., pplacer) %>%
    mutate(Species = ifelse(Relationship %in% 'Pathogenic', pplacer_sp, Species)) %>% 
    select(-id, -Relationship, -pplacer_sp) %>%
    data.frame(row.names = 1:nrow(features)) -> tax

  #  identical(rownames(dat), rownames(tax))
  
  # and parse
  ps = phyloseq(otu_table(dat, taxa_are_rows = TRUE), 
                tax_table(as(tax, 'matrix'))) 
  
  microbiome::aggregate_taxa(ps, level = 'Species')

}

prepare_ps(fileNames[1])

# consistencia de pplacer %% rdp??
features %>%
  separate(pplacer_sp, into = c('pplacer', 'prefix'), sep = ' ') %>%
  mutate(
    type = case_when(
      pplacer == Genus ~ 'Consistent',
      TRUE ~ "other"
    )) %>% 
  select(Relationship, type, Family:prefix, - Species) 

features %>% filter(grepl('Psychrobacte', pplacer_sp)) %>% select_at(ranks)



library(microbiome)
tab <- alpha(ps, index = "all")

# fantaxtic
library(fantaxtic)
# Necesitamos obtener las taxa más abundantes, en este caso el top 15
top15 <- get_top_taxa(physeq_obj = ps, n = 15, relative = T,
                      discard_other = T, other_label = "Other")
# Ya que no todas las taxa fueron clasificadas a nivel de especie, generamos etiquetas compuestas de distintos rangos taxonómicos para el gráfico
top15 <- name_taxa(top15, label = "", species = F, other_label = "Other")
# Finalmente graficamos
fantaxtic_bar(top15, color_by = "Family", 
              label_by = "Genus", facet_by = NULL, grid_by = NULL, 
              other_color = "Grey", palette = ) -> ptop15


#

# Ranked abundance distribution models for a random plot. The best model has the lowest aic.
mod <- radfit(t(ab))
mod
plot(mod, pch=".")
# or
mod <- rad.lognormal(t(ab)[1,]);plot(mod)
mod <- radfit(t(ab)[1,])
## Standard plot overlaid for all models
## Preemption model is a line
plot(mod)
## log for both axes: Zipf model is a line
plot(mod, log = "xy")
## Lattice graphics separately for each model
radlattice(mod)

# Figure: Renyi diversities in six randomly selected
# plots. The plot uses Trellis graphics with a separate
# panel for each site. The dots show the values for
# sites, and the lines the extremes and median in the
# data set.

df <- t(ab)
mod <- renyi(df)
plot(mod)
mod <- renyiaccum(df)
plot(mod, as.table=TRUE, col = c(1, 2, 2))
persp(mod)

# rarefac

rs <- rowSums(df)
quantile(rs)

Srar <- rarefy(df, min(rs))
head(Srar)
rarecurve(df, sample = min(rs))

# install.packages('iNEXT')
library(iNEXT)
# https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.html

# test <- t(ab)[1,]
# i.zero <- which(test == 0)
# test.no.zero <- test[-i.zero]

x <- apply(ab, 2, function(x) x[-which(x == 0)])

out <- iNEXT(x, q=c(0), datatype="abundance")


# Sample-size-based R/E curves, separating by "site""
ggiNEXT(out, type=1, facet.var="none", grey = T)
ggiNEXT(out, type=2, facet.var="site", grey = T)
ggiNEXT(out, type=2, facet.var="none", color.var="site")

