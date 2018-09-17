# Analisis de los metodos de clustering

El gran numero de lecturas producidas (segun la profunidad de la secuenciacion) en una corrida ofrecen una profundidad de la muestra sin precedentes. Lo que lleva a que la cola de la distribucion de abundancia de especies (biosfera 'rara') se sustancialmente mas grande y diversa de lo que anteriormente se apreciaba (sogin et al 2006 y visto en: victor kunin et al 2010, EnvMicr)

Para minimizar el error en la apreciasion de espcies (taxones) variamos el umbral de los metodos de custering debido a que no solo tomamos encuenta la varaicion evolutiva entre taxones (especies u otro nivel taxonomico defindo) si no minimizamos el efecto de la secuencacion por errores de PCR y secunciacion (Ramiro logares tal 2014, current biology).

Una de las etapas mas cruciales es el procesamiento de las lecturas para formar grupos de secuencias en base a su similitud. Estas son llamadas unidades operacionales taxonomicas (OTUs).



# Sobre-estimando filotipos raros

Cada OTU (amplicon) es interpretado com oun identificador unico de un miembro de la comunidad. 

Los errores en la secuenciacion podrian inflar la estimacion en la diversidad de la comunidad. Para tratar bioinformaticamente los errores de la secuenciacion, diversos autores han  homogenizado los siguientes parámetros (P scools, huse etal, 2007) para minimizar (asegurar) la tasa de error base-a-base (per-base) convencional de la secuenciacion Sanger, reteniendo el 90 % de los reads: (INGRESAR MEJOR TABLITA, QUE INCLUYA TIPO DE ERROR, HERRAMIENTA QUE SOLUCIONA ESTE ERROR, AUTOR)

- remover bases ambiguas (sin resolver - NNNNs) 
- Secuencias cortas o largas (aquellas secuencias en las colas de la campana de la distribucion) resolviendo de este modo errores de los primers  y barcodes.
- remover secuencias quimericas (Edgar et al)

_________



Diversos metodos han sido propuestos para agrupar (a partir de ahora usaremos la terminologia clusters o clustering) secuencias en sus 

 A continuacion, presentaremos algunas pruebas importantes desarrolladas para establecer los parametros optimos del *clustering*.

```bash
for i in $(seq 0.0 0.002 0.050);
do
echo "mothur '#system(mkdir cuttoff_$i);set.dir(input=../, output=./cuttoff_$i);set.current(current=current_files.summary);cluster.split(fasta=current, count=current, taxonomy=current, splitmethod=classify, taxlevel=7, method=opti, cutoff=$i, processors=$SLURM_NPROCS);get.current();make.shared(list=current, count=current, label=$i);classify.otu(list=current, taxonomy=current, count=current, label=$i, threshold=80)'"
done | sh
```



Y elaboramos una matriz de los resultados

```bash
for i in $(seq 0.01 0.002 0.050); do cut -f2 otus_$i/*.unique_list.0.*.cons.taxonomy | sort | uniq -c | sort -n -k2,2 | tail -n +2 | awk '{print $1}' > otus_$i.log; done
#then
paste *log > ktones.size && rm *log
sed -i 's/\t/ /g' ktones.size 
```

Y finalmente visualizamos

```bash

x <- read.csv("ktones.size", sep=" ", header=FALSE)

x[is.na(x)] <- 0

colnames(x) <- paste("CF", seq(0.01, 0.05, by = 0.002), sep="_")

pca <- prcomp(t(x), scale=FALSE)

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
 
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

pca.data <- data.frame(sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
library(ggplot2)

ggplot(data=pca.data, aes(x=X, y=Y, label=sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("Ktone-dist PCA")
```



\# okay, tienes una curva que de ascender deciente a los otus de mayor peso, me interesa que saquemos un contrapeso para visualizar 

el aumento de la pendiente usando el tamano de los otus en vez de el tamano del ktone, de este modo se remarca ambos, el peso del ktone y el tamano de la muestra (numero de otus de N ktone tamano)

O tambien podriamos ver la prueba de SAD para esta matriz de datos.

 # # Resultados

| Cutoff | # OTUs (opticlust) | % singletones |
| ------ | ------------------ | ------------- |
| 0.01   | 33862              | 63.11         |
| 0.012  | 33274              | 62.99         |
| 0.014  | 27718              | 59.63         |
| 0.016  | 27353              | 59.31         |
| 0.018  | 27023              | 59.09         |
| 0.02   | 23465              | 55.44         |
| 0.022  | 21503              | 54.04         |
| 0.024  | 21257              | 53.79         |
| 0.026  | 19925              | 52.51         |
| 0.028  | 17158              | 52.01         |
| 0.03   | 16936              | 51.74         |
| 0.032  | 15890              | 51.12         |
| 0.034  | 13540              | 50.03         |
| 0.036  | 13298              | 49.65         |
| 0.038  | 12970              | 49.28         |
| 0.04   | 11669              | 48.31         |
| 0.042  | 11249              | 48.20         |
| 0.044  | 10994              | 47.76         |
| 0.046  | 10357              | 47.29         |
| 0.048  | 9512               | 47.51         |
| 0.05   | 9341               | 47.65         |

 