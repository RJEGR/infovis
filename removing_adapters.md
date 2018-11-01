El objetivo de esta prueba es ver el efecto que tiene remover o conservar los primers en ambos, la base de referencia, asi como las lecturas de cada una de las bibliotecas secuencidas. A continuacion se muestra el el archivo fasta de ambos adaptores+ primer en su direccion forward y reverse (se subraya unicamente el primer)

>AdV_1389F length=48
>TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG**TTGTACACACCGCCC**
>AdV_1510R length=54
>GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG**CCTTCYGCAGGTTCACCTAC**

Podemos visualizar el numero de lecturas que potencialmente conservan el primer:

```bash
zcat 4-X04-A10-18S-AMB_S7_L001_R1_001.fastq.gz | head | grep "TTGTACACACCGCCC" --color
```

En el siguiente ejemplo encontramos un Total seq = 138561 , de las cuales seq con primer 135483 <- el resto de las seq son de mala calidad al parece

The Forward primer looks like:

```bash
for i in $(ls R1*fastq.gz);
do
echo $i $(zcat $i | grep -c "TTGTACACACCGCCC") $(zcat $i | grep -c "^@");
done

# second paired-end
for i in $(ls R2*fastq.gz);
do
echo $i $(zcat $i | grep -c "CCTTC[A,C,T,G]GCAGGTTCACCTAC") $(zcat $i | grep -c "^@");
done
```

## Cutadapter

Using cutadapter to remove partial adapters. If you have paired-end data, trim also read 2 with the revere complement of the Truseq Universal Adapter. The full command-line will run as follow:

```bash
srun ./cutadapt \
            -j 24 \
            -g TTGTACACACCGCCC \
            -G CCTTCYGCAGGTTCACCTAC \
            -o trimmed.1.fastq.gz -p trimmed.2.fastq.gz \
            04-X04-C21-18S-AMB_S17_L001_R1_001.fastq.gz 04-X04-C21-18S-AMB_S17_L001_R2_001.fastq.gz &

```

> Cutadapter RECONOCE PARCIALMENTE ADAPTORE DEBIDO A QUE RECONOCE EL FORMATO Y, ETC., [Ver manual](https://cutadapt.readthedocs.io/en/stable/guide.html#read-processing-stages) para mayor detalles.

En pruebas anteriores se habia utilizado una base de referencia curada (reference=silva132_kp_aln.fasta,) en la que se incluian la region que flanquean los primers utilizando los siguientes parametros: screen.seqs(fasta=current, count=current, summary=current, minlength=96, start=17, end=1426, maxhomop=10).

>  resultados en /LUSTRE/bioinformatica_data/genomica_funcional/mothur/XIXIM0_old) 

La siguiente prueba nos arroja los siguientes parametros a utilizar despues de utilizar las lecturas ore-recordadas, generar los contigs, y alinearla con la el modelo de ssu y la base de referencia silva132_ep_v9_aln.fasta:

```bash
screen.seqs(fasta=current, count=current, summary=current, minlength=90, start=17, end=1391, maxhomop=10)
```

Filtramos datos en bash

```bash

awk '{print $4, "Contigs"}' ./CONTIGS/cigom.trim.contigs.summary > cigom.trim.contigs.summary.lenghts
awk '{print $4, "GContigs"}' ./CLEAN_CONTIGS/cigom.trim.contigs.good.summary > cigom.trim.contigs.good.summary.lenghts
awk '{print $4, "Unique"}' ./CLEAN_CONTIGS/cigom.trim.contigs.good.unique.pick.summary > cigom.trim.contigs.good.unique.pick.summary.lenghts
awk '{print $4, "SSU"}' ./ALIGNED_CONTIGS/cigom.trim.contigs.good.unique.pick.summary > cigom.trim.contigs.good.unique.pick.summary2.lenghts
awk '{print $12, "Aligned"}' ./ALIGNED_CONTIGS/cigom.trim.contigs.good.unique.pick.align.report > align.report.lenghts

```

Y renombramos manualmente las cabeceras de los archivos utilizando el titulo "nbases" y "Step" para cada columna.

Cargamos en R

```R
tmp <- list.files(path = "./", pattern = "*.lenghts", full.names = TRUE)
myfiles = lapply(tmp, read.csv, sep = " ", header=TRUE, stringsAsFactors=FALSE)
data <- data.frame(do.call(rbind, myfiles))

library(ggplot2)
#ggplot(data, aes(nbases, fill = step)) + geom_density(alpha = 0.2, aes(y = ..count..))

p <- ggplot(data, aes(nbases, fill = Step)) + geom_histogram(alpha = 0.7, aes(y = ..ncount..), position = 'identity') + theme_classic()
p + scale_fill_brewer(direction = -1, palette = "Set1")
 
            #facet_wrap(~ Step) + 
            #scale_y_continuous(trans = 'log10')

#p <- ggplot(data, aes(nbases, colour = step)) + geom_freqpoly(alpha = 1, aes(y = ..count..), position = 'identity') + theme_classic()
#p + scale_color_brewer(direction = -1, palette = "Set1")

# :::
hist(myfiles[[1]][,1], col='skyblue',border=F)
hist(myfiles[[2]][,1],add=T,col=scales::alpha('red',.5),border=F)
# :::
library(psych)
py <- data.frame(align=myfiles[[3]][,1], unique=myfiles[[2]][,1])
multi.hist(py, dcol = "red")

# :::
pairs.panels(py, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
             )
```

> please, refers to the nanodrop plot example showed [here](https://gigabaseorgigabyte.wordpress.com/2017/06/01/example-gallery-of-nanoplot/) to the final pairs.panel plot



Veamos que hay con esas secuencias que estan por debajo del threshold de tamano. Usamos el ejemplo de filtrado con awk visto [aquí](https://www.tim-dennis.com/data/tech/2016/08/09/using-awk-filter-rows.html)

```bash
awk '{if ($4 <= 90) {print} }' cigom.trim.contigs.good.unique.pick.summary > cigom.trim.contigs.good.unique.pick.summary.90less
awk '{print $4}' cigom.trim.contigs.good.unique.pick.summary.90less | sort | uniq -c | sort -k2,2      
	  1 14
      1 22
      7 3
      2 37
      1 4
      4 5
      1 56
      1 6
      1 8
      2 86
      6 87
      2 88
      4 89
      1 90
```

## Analizando los resultados del alineamiento



```R
x <- read.csv("cigom.trim.contigs.good.unique.pick.align.report", sep="\t")
str(x[1,])
'data.frame':   1 obs. of  16 variables:
 $ QueryName              : Factor w/ 132043 levels "M03978_33_000000000..."
 $ QueryLength            : int 137
 $ TemplateName           : Factor w/ 1463 levels "AACY020379403.855.2641"
 $ TemplateLength         : int 130
 $ SearchMethod           : Factor w/ 1 level "kmer": 1
 $ SearchScore            : num 94.6
 $ AlignmentMethod        : Factor w/ 1 level "needleman": 1
 $ QueryStart             : int 2
 $ QueryEnd               : int 131
 $ TemplateStart          : int 1
 $ TemplateEnd            : int 130
 $ PairwiseAlignmentLength: int 130
 $ GapsInQuery            : int 0
 $ GapsInTemplate         : int 0
 $ LongestInsert          : int 0
 $ SimBtwnQuery.Template  : num 100
```



En este ejemplo, vemos que la secuencia candidato _M03978_33_000000000..._  tenía una longitud de **137** **bases** y usamos la búsqueda **kmer** (es decir, usando 8mers) para encontrar la secuencia de plantillas _AACY020379403.855.2641_ como la mejor coincidencia. La secuencia _AACY020379403.855.2641_ tuvo el **94.6%** de los 8mers encontrados en la secuencia candidato.

A continuación, vemos que usamos el método de alineación de **_Needleman_** resultó en una alineación de pares de 131 caracteres. 

En este ejemplo, vemos que la alineación realmente comienza en la posición de secuencia candidata 2. Esto ocurrió porque la secuencia realmente tiene una secuencia vectorial en el extremo 5 ', que no está representada en la alineación de referencia, que incluye los sitios de cebadores tradicionales 27f y 1492r. 

A continuación, vemos que durante la etapa de alineación por pares, se introdujeron 0 espacios tanto en la secuencia candidata como en la secuencia de la plantilla; ninguno de los espacios de la plantilla tuvo que ser corregido con el algoritmo NAST. Finalmente, la secuencia candidata alineada fue idéntica en un 100% a la secuencia de la plantilla (incluidos los huecos) ([mothur 2018](https://www.mothur.org/wiki/Align.seqs))

> Correlacionando datos, 

```R
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
M <- cor(x[,-c(1,3,5,7)], method = "pearson")

p.mat <- cor.mtest(x[,-c(1,3,5,7)])
library(corrplot)
# install.package(corrplot)
corrplot(M, type="upper", order="hclust", col=c("black", "white"),
         bg="lightblue",
         tl.col="black", tl.srt=45,
         p.mat = p.mat, sig.level = 0.01, insig = "blank")
```



Y en base. a estos resultados, podemos revisar algunas relaciones mas estrechas

```R
library(psych) # high demand time
pairs.panels(x[, c(2,6,12,16)], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
             )

```

Y encontramos que a medida que en tamano del alineamiento es, el % de similitud es mayor

O replicar la visualacion de [aqui]()

```R
dev.off()
# I divide the screen in 2 line and 1 column only
my_screen_step1 <- split.screen(c(2, 1))
# I add one graph on the screen number 1 which is on top :
screen(my_screen_step1[1])
plot( x[,12],x[,16] , pch=20 , xlab="Tamano del alineamiento" , ylab="% Similitud", cex=3 , col=rgb(0.4,0.9,0.8,0.5) )
 
# I divide the second screen in 2 columns :
my_screen_step2=split.screen(c(1, 2), screen = my_screen_step1[2])
screen(my_screen_step2[1])

hist(x[,12], border=F , col=rgb(0.2,0.2,0.8,0.7) , main="" , xlab="Distribucion del Tamano del alineamiento")
screen(my_screen_step2[2])
hist(x[,16], border=F , col=rgb(0.8,0.2,0.8,0.7) , main="" ,  xlab="% De similitud")

```

Un paso adicional

```R

# data <- data.frame(x = x$PairwiseAlignmentLength, y = x$SimBtwnQuery.Template)
data <- data.frame(x = x[,2], y = x[,16])

d <- ggplot(data, aes(x=x, y=y) )
d + stat_density_2d(geom = "raster", aes(fill = ..density..), contour = FALSE)

# OR
m <-  ggplot(data, aes(x=x, y=y) ) + geom_point()
m + stat_density_2d(aes(fill = stat(level)), geom = "polygon")
```

O

```r
ggplot(x, aes(PairwiseAlignmentLength, TemplateLength)) + geom_point() + geom_rug(col = "darkred", alpha = 0.1)

```

O mejor

```bash

data <- data.frame(x = x[,2], y = x[,16])

# placeholder plot - prints nothing at all
empty <- ggplot() + geom_point(aes(1, 1), colour = "white") + theme(plot.background = element_blank(),
panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.border = element_blank(), panel.background = element_blank(), axis.title.x = element_blank(),
axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),
axis.ticks = element_blank())

#Assign color variables
col1 = "#d8e1cf" 
col2 = "#438484"

# Density of x and y variables
p1 <- ggplot(data, aes(x=x,y=y))+
  stat_density2d(aes(fill=..level..), geom="polygon") + 
  scale_fill_gradient(low = col1, high = col2)  + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        #axis.ticks = element_blank(),
        panel.background = element_blank()) +
        coord_cartesian(xlim = c(130, 150), ylim = c(75, 100)) # para la EP
        # coord_cartesian(xlim = c(160, 175), ylim = c(80, 100)) <-- para KP 


p2 <- ggplot(data, aes(x = x)) + stat_density(fill = col1) + geom_rug(col = col1, alpha = 0.1) + theme(panel.background = element_blank()) + xlab("Tamano de lectura (query)") + coord_cartesian(xlim = c(130, 150))
p3 <- ggplot(data, aes(x = y)) + stat_density(fill = col2) + geom_rug(col = col2, alpha = 0.1) + coord_flip(xlim = c(75, 100)) + theme(panel.background = element_blank()) + xlab("Identidad del Alineamiento") 

# arrange the plots together, with appropriate height and width for each row
# and column

library(gridExtra)
grid.arrange(p2, empty, p1, p3, 
                 ncol=2, 
                 nrow=2, 
                 widths=c(4, 1), 
                 heights=c(1, 4))


```

Segundo visualizacion

```R
data <- data.frame(x = x[,4], y = x[,12], Identidad = x[,16])

#Assign color variables
col1 = "#d8e1cf" # "#0091ff"
col2 = "#438484" # "#f0650e"

# 
library(ggpubr)
p <- ggplot(data, aes(x=x, y=y, color=Identidad)) +
      geom_point() +
      scale_color_gradient(low = col1, high = col2) + scale_alpha(range = c(.05, .25)) +
      theme(legend.position="bottom",
            panel.background = element_blank()) +
      xlab("Tamano de lectura (query)") + ylab("Tamano del alineamiento") + 
      stat_cor(method = "pearson")

# devtools::install_github("daattali/ggExtra") 
library(ggExtra)

# with marginal histogram
ggMarginal(p, type="histogram")
# marginal density
ggMarginal(p, type="density")
# marginal boxplot
ggMarginal(p, type="boxplot")


```

# Comprobar adaptores en las bases de referencia

En el laboratorio utilizamos dos bases de datos para alinear nuestros contigs, estan son la silva 132 en las que se incluye y no incluye los fragmentos de los primers (base kp y ep, respectivamente). Para corroborar que los primers se conservan o han sido removidos en ambas bases de referencia, utilzaremos cutadapt del mismo modo que lo hicimos con las bibliotecas crudas.

Primero preparamos las bases que estan alineadas del siguiente modo:

```bash
cat silva132_kp_aln.fasta | sed 's/-//g' |awk '{ print toupper($0) }' > silva132_kp.fasta

```

Y ejecutamos cutadapter con la base kp

```bash

srun ./cutadapt -g TTGTACACACCGCCC -a CCTTCYGCAGGTTCACCTAC -o silva132_kp.fasta.cutadapt silva132_kp.fasta 2> silva132_kp_aln.fasta.cutadapt.log &

# siguiente paso
cat silva132_kp.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > silva132_kp.fasta.lengths 

# siguiente paso
cat silva132_kp.fasta.cutadapt | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > silva132_kp.fasta.cutadapt.lengths

# finalmente
paste *lengths | awk '{print $1, $2, $4}' > kp.lengths
```

> Regular 3’ adapter	-a ADAPTER
> Regular 5’ adapter	-g ADAPTER

Repetimos lo mismo para la siguiente base (ep):

```bash
cat silva132_ep_v9_aln.fasta | sed 's/-//g' |awk '{ print toupper($0) }' > silva132_ep_v9.fasta

srun ./cutadapt -g TTGTACACACCGCCC -a CCTTCYGCAGGTTCACCTAC -o silva132_ep_v9.fasta.cutadapt silva132_ep_v9.fasta 2> silva132_ep_v9.fasta.cutadapt.log &

# siguiente paso:
cat silva132_ep_v9.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'> silva132_ep_v9.fasta.lengths

# siguiente paso:
cat silva132_ep_v9.fasta.cutadapt | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > silva132_ep_v9.fasta.cutadapt.lengths

# finalmente
paste silva132_ep_v9.fasta.cutadapt.lengths silva132_ep_v9.fasta.lengths | awk '{print $1, $2, $4}' > ep.lengths

```

Y cargamos en R

```R
kp <- read.csv("kp.lengths", sep=" ", header=TRUE)
colnames(kp) <- c("id", "cutadapt", "raw")

ep <- read.csv("ep.lengths", sep=" ", header=TRUE)
colnames(ep) <- c("id", "cutadapt", "raw")



```

Y Visualizamos su distribucion 

```R
ep$Reference <- "EP"
kp$Reference <- "KP"
dataplot <- rbind(ep, kp)
dataplot <- melt(dataplot)

p <- ggplot(dataplot, aes(value, fill = variable)) + geom_histogram(alpha = 0.7, bins = 200, aes(y = ..ncount..), position = 'identity') + theme_classic()
p + scale_fill_brewer(direction = -1, palette = "Set1") + facet_wrap(~ Reference)
```

![](/Users/cigom/Desktop/Screen Shot 2018-10-24 at 2.39.20 PM.png)