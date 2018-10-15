# Introduccion

El siguiente documento describe las pruebas bionformaticos realizados durante el procesamiento de amplicones 18S - region V9. 

Los resultados descritos justifican cada uno de los parametros definidos, de este modo se transfiere al usuario la mejor practica en el manejo de amplicones y se deja al alcance una serie de flujos de trabajo para replicar las pruebas bioinformaticas o como punto de referencia para analisis metagenomicos similares.

En este documento se obvia la implementacion de bibliotecas paired-end generadas por un secuenciador tipo Illumina MiSeq.  La tubería del análisis implementada en esta sesion es modificada del sugerido por Schloss et al. 2013 para amplicones de la región v9 del gen 18S rRNA secuenciados por Illumina obtenidas de alicuotas de origen marino.

disculpen la falta de acentos debido a la incompatiblidad de caracteres especiales de mi teclado.

- Bases de referencia
  - curacion y enriquecimiento de bases de referencia
- procesamiento de las lecturas
  - generando contigs
  - Colapsando secuencias repetidas
  - alineando a base de referencia
  - ajustando al modelo probabilistico de covariancia de la region ribosomal de interes
- remover secuencias quimericas
  - agrupando contigs para remover quimeras
- Asignacion taxonomica
  - bases de referencia curadas
    - la importancia de un marco de referencia taxonomico
  - problemas y retos https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3501-4
- Generando unidades operacionales Taxonomicas
  - segundo agrupamiento de las secuencias 
  - reasignacion taxonomica de OTUs
  - radio de inflacion taxonomica (LULU)
- Analisis cascada abajo
  - generar base de datos de resultados  y otros archivos para analisis de cascada abajo
  - Procesamiento de taxones raros (1-tones)
  - Inferencias estadisticos en base a hipotesis nulas
  - Correlacion de taxones con alguna variable ambiental
  - Visualization de resultados

# Preparando los archivos para analisis cascada abajo

- phyloseq
- base de datos

## Creando una base de datos del analisis metagenomico

Just sending this email to note a bug when trying to create database from my metagenomic analysis.

According to the description of the create.database ([default settings](https://mothur.org/wiki/Create.database)) if you set other parameter configuration in the get.oturep, the create.database module does not work properly, printing an error as a follow example:

OTU size info does not match for bin 3. The contaxonomy file indicated the OTU represented 37 sequences, but the repfasta file had 1.  These should match. Make sure you are using files for the same distances

The configuration from my command than print in this error is as follow:

```bash

system(mkdir database)
set.dir(input=../MAKE_OTUS, output=./database)
set.current(current=current_files.summary)

get.oturep(list=current, label=$CUTOFF, fasta=current, name=current, method=abundance)
classify.otu(list=current, name=current, taxonomy=current, label=$CUTOFF)

create.database(list=current, label=$CUTOFF, repfasta=./database/${PROJECT}.trim.contigs.good.unique.pick.good.filter.unique.precluster.pick.opti_mcc.unique_list.${CUTOFF}.rep.fasta, repname=./database/${PROJECT}.trim.contigs.good.unique.pick.good.filter.unique.precluster.pick.opti_mcc.unique_list.${CUTOFF}.rep.names, constaxonomy=./database/${PROJECT}.trim.contigs.good.unique.pick.good.filter.unique.precluster.pick.opti_mcc.unique_list.{$CUTOFF}.cons.taxonomy)


```





Finaly replying the default settings in the **get.oturep**(list=current, label=0.03, fasta=current, name=current) my database was created properly !. I hypothesize than the error here could be due to the sorted way of the fasta file when switch on the option sorted=size in get.oturep, Do you think is necessary mention it or fix it in the [manua](https://mothur.org/wiki/Create.database)l?

By the way, could you include in your newer version of mothur the repfasta and repname or repcount in the current files in order to save the current names in the set.current command? 

Regards!

## Cribando secuencias

While many rare OTUs are beyond doubt real biological entities[10](https://www.nature.com/articles/s41467-017-01312-x#ref-CR10), an appreciable fraction of rare OTUs are likely errors from PCR and sequencing[11](https://www.nature.com/articles/s41467-017-01312-x#ref-CR11), (lulu paper)

Before estimating the number of OTUs predicted by pyrosequencing of the pooled template preparations, we removed low‐quality reads and identified sequences that did not represent targeted hypervariable regions. In the case of DNA extracted from *E. coli* and *S. epidermidis*cultures, **many of the non‐target reads mapped to specific regions of the genome that did not code for rRNAs. In some genomes our rRNA primers will bind with low efficiency to similar but non‐identical targets to produce non‐rRNA pyrotags.**

After removing low‐quality sequences (Huse *et al.*, 2007; Kunin *et al.*, 2010) and reads that did not represent the targeted region, we measured the per‐base error rate to be 0.0021–0.0042. ([ref](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1462-2920.2010.02193.x))



# tratar singletones (artefactos?)

Our data indicate that many singletons (average of 38 % across gene regions) are likely artifacts or potential artifacts, but a large fraction can be assigned to lower taxonomic levels with very high bootstrap support(w32 % of sequences to genus with 90 % bootstrap cutoff). Further, many singletons clustered into rare OTUs from other datasets highlighting their overlap across datasets or the poor performance of clustering algorithms. 

 

OTUs were considered artifacts if: 

1. OTUs were unclassiﬁed at a phylum level (many uncultured sequences may lack phylum level classiﬁcation thus exaggerating proportion of artifact OTUs); 
2. they did not classify to a phylum at 50 % bootstrap support or higher; or, 
3. the ITS sequences could not be mapped to ITS1 or ITS2 region (ITSx)

 To investigate if these singletons represent true biological or artiﬁcial variability (platform speciﬁc variability, indels due to polymerase slippage, or homopolymeric reads), we aligned singletons against representative sequences of the 100 most abundant OTUs from the original datasets. The mismatches among singletons and the representative sequences of the common OTUs generated on 454 and Illumina platforms appeared stochastically distributed across the alignments suggesting that they were unlikely a result of poor read quality in the read termini. Singletons generated using 454 tech- nologies differed from abundant OTUs frequently because of inconsistent homopolymer lengths and/or single nucleotide differences. In contrast to 454-sequencing, differences in the Illumina-generated singletons were most often nucleotide differences with no evidence of inconsistent homopolymer lengths. Based on these ﬁndings it is impossible to determine the source of the variability as polymerase slippage, sub-optimal platform performance or true biological variability could result in similar outcomes

 

(PDF) Scraping the bottom of the barrel: Are rare high throughput sequences artifacts?

https://www.researchgate.net/publication/266562420_Scraping_the_bottom_of_the_barrel_Are_rare_high_throughput_sequences_artifacts#pf5

**Scraping the bottom of the barrel: Are rare high throughput sequences artifacts?**

https://www.researchgate.net/publication/266562420_Scraping_the_bottom_of_the_barrel_Are_rare_high_throughput_sequences_artifacts



# Segundo agrupamiento de las secuencias 

Escribir resumen del siguiente articulo. https://aem.asm.org/content/79/21/6593 ref figura 1 sobre distribucion a lo largo de muestras.

## LULU: Metodo de pos-agrupamiento de amplicones basado en la distribucion de OTUs

The algorithm is intended as a post-clustering OTU table curation method aimed at removing erroneous OTUs from tables produced by any clustering algorithm e.g., methods used in this study[13](https://www.nature.com/articles/s41467-017-01312-x#ref-CR13),[22](https://www.nature.com/articles/s41467-017-01312-x#ref-CR22),[23](https://www.nature.com/articles/s41467-017-01312-x#ref-CR23),[24](https://www.nature.com/articles/s41467-017-01312-x#ref-CR24), and those implemented in Qiime[28](https://www.nature.com/articles/s41467-017-01312-x#ref-CR28) and Mothur[29](https://www.nature.com/articles/s41467-017-01312-x#ref-CR29), as long as the product is an OTU table and a corresponding file with representative sequences.

The implementation of the algorithm is based on a set of assumptions based on four observations we have previously made when working with HTS of amplified marker genes (a.k.a. metabarcoding) of well-studied organism groups with well-populated reference databases present (i.e., plants, as used for validation here). The first observation is that OTU tables often have more OTUs than expected from biological knowledge of the system under investigation[11](https://www.nature.com/articles/s41467-017-01312-x#ref-CR11). The second observation is that OTU tables often contain low-abundance OTUs, which are taxonomically redundant in the sense that their taxonomic assignment is identical to more abundant OTUs. This pattern may be caused by incomplete reference data and/or insufficient clustering, but can also indicate that the OTU is effectively a methodological artefact. The third observation is that the highest sequence similarity (match rate) of such taxonomically redundant, low-abundance OTUs with any reference sequence is most often low compared to the sequence similarity of more abundant OTUs with the same taxonomic assignment. The fourth observation is that such seemingly redundant and less abundant OTUs almost consistently co-occur (i.e., are present in the same samples) with more abundant OTUs with a better taxonomic assignment. Based on these observations, it can be assumed that the majority of these low-abundant OTUs are in fact methodological and/or analytical errors, or rare (intragenomic) variants, which will cause inflated diversity metrics. Following from this assumption, the LULU algorithm is constructed to iteratively work though the OTU table to flag potential erroneous OTUs by employing the observed patterns of co-occurrence guided by pairwise similarity of centroid sequences of the OTUs. Thus, the algorithm takes advantage of the observed reproducible nature of extra/spurious OTUs and their sequence similarity to more abundant OTUs in the same samples and uses these features to infer their nature as errors (or true—but taxonomically redundant) variants of biological entities already represented in the table. After identification of these extra OTUs, they can be merged with their parent OTUs in order to preserve the total read count and reduce the OTU number of the table to a biologically reasonable level. The resulting table may be subjected to direct species richness metrics and other biodiversity analyses dependent on species-level OTU delimitation.

```bash
# get otu rep for LULU TEST BASED ON DISTANCE
for i in $(seq 0.010 0.002 0.050); 
do 
echo "mothur '#system(mkdir lulu_$i);set.dir(input=./otus_$i, output=./lulu_$i);set.current(current=current_files.summary);get.oturep(column=current, label=$i, fasta=current, name=current, method=distance)'"
done | sh


```

The match list can in practice be produced with any tool for pair wise matching of sequences. **BLASTn** is an effective tool for this, and the one that was used in the validation of the **LULU** algorithm. The only requirements is that the match list has three columns with pair wise similarity scores for the OTUs. The first column contains the id of the query OTU, the second column contains the id of the matching OTU, and the third column contains the similarity score (%) of the two OTUs.

The match list can be produced with **VSEARCH** with these commands:

 ```bash
# 1) Parsing the representative sequence with otu labels
# 2) removes columns from alignments based on a criteria defined by every character is either a '.' or a '-' 

for i in $(seq 0.010 0.002 0.050); 
do
awk '/^>/{gsub(/[|]/, " "); print ">"$2; next}{print}' lulu_$i/cigom.trim.contigs.good.unique.pick.good.filter.unique.precluster.pick.opti_mcc.unique_list.*.rep.fasta | sed 's/[-, .]//g' > lulu_$i/dist.${i}.rep.fasta;
done

# 3) Then, produce the match list

for i in $(seq 0.010 0.002 0.050); do srun vsearch --usearch_global lulu_$i/dist.${i}.rep.fasta --db lulu_$i/dist.${i}.rep.fasta --self --id .84 --iddef 1 --userout vsearch.${i}.list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10 --threads 96; done | sh 2> vsearch.log &
# 4) 
mkdir vsearch
mv vsearch.*.list.txt ./vsearch
 ```

BLASTn - Alternatively, a matchlist can also be produced with **BLASTn** with this command  (slow process!!). Ex.

```bash
makeblastdb -in lulu_$i/dist.0.010.rep.fasta -parse_seqids -dbtype nucl

blastn -db lulu_$i/dist.0.010.rep.fasta -outfmt '6 qseqid sseqid pident' -out blastn_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query lulu_$i/dist.0.010.rep.fasta

```



Esta semana (y parte de la anterior) estuve probando algunas cosas:

1. Lulu: Cribado de OTUs (potencialmente erróneos) a través de las matrices de co-ocurrencia a lo largo de un grupo de muestras (sea nuestro caso las estaciones de c/ crucero) y matrices de distancia de la similitud entre secuencias (usando vsearch y blastn, o la estrategia de upticlust). 

2. Selección de secuencia representativa de dos maneras:

   1. Distancia de similitud media mínima del cluster de secuencias (Por OTU) ← 

   2. Secuencia representativa por su abundancia (ie. más abundante)Se seleccionó la secuencia representativa del punto (1) debido a que más tarde son comparadas las secuencias de c/ uno de los OTUs entre sí.Esto nos llevó a dos cosas:

      - Implementar vsearch / blastn para re-calcular la similitud entre secuencias a través de un alineamiento de pares de bases (la estrategia de vsearch fue la implementada debido a la rapidez con la que esta concluye)

      - Recuperar la distancia previamente calculada durante el análisis con upticlust para cada una de esas secuencias representativas usando el módulo de mothur get.dist.

      - Comparación de test de mantel par determinar variación de las matrices

        

        Visualización de grafos de distancia de la similitud de secuencias (matrices cuadráticas de distancia)

        

      - Debido a que la serie de datos de distancia es muy grande buscamos reducir las dimensiones basado en la hipótesis de que los datos se agruparán de acuerdo a clusters de los OTUs.  Método PCA y tSNEA pesar de reducir las dimensiones. 

      - La visualización de tantas datos es muy tardado, entonces descifremos el modo de reducir la cantidad de información.El objetivo de lulu es cribar con aquellos taxones raros que inflan la estimación de la alfa diversidad, y su distribución suele estar en alas colas, es decir los k-tones, especialmente los singletones. 

      - La hipótesis es que los singletones se reducirán al “reintegrarse” a los otros abundances en el procesamiento con Lulu (ver funcionamiento del algoritmo). Resultados:Usando el segundo agrupamiento con vsearch se logró visualizar los grafos coloreando los nodos por su clasificación taxonómica. Y haciendo un facet_wrap cada Reino. Esto demora 15-20 minutos en procesarse con aprox. 200 K nodos. 

      - Sin embargo no se revela algun patron interesante.Se visualizó resultados de la reducción de dimensiones con tSNE usando un subset de tres grupos de OTUs (resultado vsearch) pero sin algun patron satisfactorio.Usando 3 OTUs (nuevamente, resultados vsearch). Se logró visualizar el grafo de similitud , en el que la información interactúa de acuerdo a su grupo (ie. OTUs al que corresponde dicho nodo (cada nodo representa una secuencia agrupada en los OTUs); también, se logró etiquetar c/ su asignación a N nivel taxonómico. los resultados muestran 3 grupos ideales de OTUs, sin embargo la clasificación es heterogénea a lo largo de los nodos que se agrupan.Se comparó la alfa diversidad y se demuestra una correlación positiva entre los datos procesados con lulu y los datos ‘crudos’ de mothur.

# La importancia de un marco de referencia taxonomica

Most life on earth is microbial, belonging to the ‘Bacteria’ and ‘Archaea’ domains (1), and to numerous lineages of microbial ‘Eukaryota’ (e.g. protists) (2). Less than 1% of microbes are cultivable, and therefore diversity was vastly underestimated by traditional microbiological methods (3). The known extent of microbial diversity has grown and continues to grow rapidly as sequence-based methods are used to characterize microbes (4). One of the major breakthroughs in the study of the diversity of microbes was the use of the ribosomal rRNA (rRNA) gene sequences, particularly of the small subunit (SSU; also called 16S rRNA for Bacteria and Archaea and 18S rRNA for Eukaryota) (https://academic.oup.com/nar/article/42/D1/D643/1061236)