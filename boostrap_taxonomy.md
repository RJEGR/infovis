# Efecto en el boostrap

Tomamos como ejemplo dos set de datos que fueron procesados con distintas bases de referencia (EP y KP). Comparemos su efecto en el boostrap.

```bash
PATH=/LUSTRE/bioinformatica_data/genomica_funcional/mothur/X04_mothur_18S_05102018_EP
# y 
/LUSTRE/bioinformatica_data/genomica_funcional/mothur/X04_mothur_18S_10092018/CLASSIFY

scp rgomez@omica:$PATH/*.taxonomy .
```

Abrimos los datos en r:

```R


kp <- read.csv("cigom.trim.contigs.good.unique.pick.good.filter.unique.precluster.pick.w2pr2_worms_API02.wang.taxonomy.KP", sep="\t", header=TRUE, stringsAsFactors=FALSE)
ep <- read.csv("cigom.trim.contigs.good.unique.pick.precluster.pick.w2pr2_worms_API02.wang.taxonomy.EP", sep="\t", header=TRUE, stringsAsFactors=FALSE)

constaxa <- kp

tax <- strsplit(constaxa[,2], ";")
tax <- sapply(tax, "[", c(1:7))
tax <- as.data.frame(t(tax))

Ranks <- c("Reino", "Filo", "Clase", "Orden", 
                   			 "Familia", "Genero", "Especie")

# ::::::::::::::::::: Remove non-numeric character by:  <- EN CASO DE  UN ANALISIS DEL BOSTRAP DE LA ASIGNACION
Q <-as.data.frame(apply(tax, 2, 
    function(x) gsub("[^0-9\\.]", "",  x, perl=TRUE)),
     stringsAsFactors = F)

Q <- apply(Q, 2, 
    function(x) as.numeric(x))
Q[is.na(Q)] <- 0
colnames(Q) <- Ranks


Qep <- Q
Qkp <- Q



library(reshape2)
v1 <- melt(Qkp, 
          id.vars = colnames(Qkp),
          value.name= c("Pct"))

v2 <- melt(Qep, 
          id.vars = colnames(Qep),
          value.name= c("Pct"))


v1$Reference <- "KP"
v2$Reference <- "EP"


violin <- rbind(v1,v2)
violin$Reference <- as.factor(violin$Reference)

library(ggpubr)
ggviolin(violin, x = "Var2", y = "Pct", 
                    fill = "Reference", color = "Reference", 
                    add = "none", add.params = list(fill = "white")) + 
                    xlab("Nivel Taxonomico") +
                    ylab("% de asignacion") 

# =========

p <- ggboxplot(violin, x = "Reference", y = "Pct",
                color = "Reference", palette =c("#00AFBB", "#E7B800"),
                add = "none", shape = "Reference", facet.by = "Var2", short.panel.labs = TRUE) +
                ylab("% de asignacion") + xlab("Nivel Taxonomico") 

p + stat_compare_means(comparisons = list(c("KP", "EP")), label.y = c(25, 50, 100))+
stat_compare_means(label.y = 45)  




ggplot(violin, aes(Pct, ..count.., colour=Var2, fill=Var2)) +
    geom_density(alpha=0.35, position = "fill") +
    scale_fill_brewer(palette = "Paired" ) +
    scale_color_brewer(palette ="Paired" ) +
    facet_wrap(~ Reference) + theme_minimal()

ggplot(violin, aes(Pct, ..count.., colour=Reference, fill=Reference)) +
    geom_density(alpha=0.35, position = "stack") +
    scale_fill_brewer(palette = "Paired" ) +
    scale_color_brewer(palette ="Paired" ) +
    facet_wrap(~ Var2) + theme_minimal()

           
ggplot(violin, aes(Pct, ..count.., colour=Reference, fill=Reference)) +
    geom_histogram(alpha=0.35, position = "stack", bins = 100) +
    scale_fill_brewer(palette = "Paired" ) +
    scale_color_brewer(palette ="Paired" ) +
    facet_wrap(~ Var2, scales="free_y") + theme_minimal()           


```

![](/Users/cigom/Desktop/Screen Shot 2018-10-24 at 1.03.47 PM.png)

Y ademas,

![](/Users/cigom/Desktop/Screen Shot 2018-10-24 at 2.47.14 PM.png)Y otro tipo

![](/Users/cigom/Desktop/Screen Shot 2018-10-24 at 1.00.49 PM.png)