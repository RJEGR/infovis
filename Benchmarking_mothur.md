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

for i in $(seq 0.010 0.002 0.050); do echo "srun vsearch --usearch_global lulu_$i/dist.${i}.rep.fasta --db lulu_$i/dist.${i}.rep.fasta --self --id .84 --iddef 1 --userout vsearch.${i}.list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10"; done | sh &> vsearch.log &

# 
mkdir vsearch
mv *list.txt ./vsearch
 ```

BLASTn - Alternatively, a matchlist can also be produced with **VSEARCH** with this command  (slow process!!)

```bash
makeblastdb -in phylo.rep.treein.fasta -parse_seqids -dbtype nucl

blastn -db phylo.rep.treein.fasta -outfmt '6 qseqid sseqid pident' -out blastn_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query phylo.rep.treein.fasta

```





# La importancia de un marco de referencia taxonomica

Most life on earth is microbial, belonging to the ‘Bacteria’ and ‘Archaea’ domains (1), and to numerous lineages of microbial ‘Eukaryota’ (e.g. protists) (2). Less than 1% of microbes are cultivable, and therefore diversity was vastly underestimated by traditional microbiological methods (3). The known extent of microbial diversity has grown and continues to grow rapidly as sequence-based methods are used to characterize microbes (4). One of the major breakthroughs in the study of the diversity of microbes was the use of the ribosomal rRNA (rRNA) gene sequences, particularly of the small subunit (SSU; also called 16S rRNA for Bacteria and Archaea and 18S rRNA for Eukaryota) (https://academic.oup.com/nar/article/42/D1/D643/1061236)