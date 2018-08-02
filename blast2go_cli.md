# Blast2go_cli - testing 

> Ricardo Go-Re

## Create a Fasta file database for local Blast and to import XML results successfully into Blast2GO

> https://www.blast2go.com/support/blog/22-blast2goblog/111-format-fasta-file-blast

The Blast2GO Command Line is a professional solution for flexible, high-performance and automatic functional annotation tasks. This Annotation Pipelineallows your to integrate and automate your functional annotation task in a flexible way (Blast2GO Command Line Use Manual, 2015). The Command Line is based on the Blast2GO methodology, first published in 2005, [Conesa et al., 2005] for the functional annotation and analysis of gene or protein sequences. The method uses local sequence alignments (BLAST) to find similar sequences (potential homologous) for one or several input sequences. The program extracts all Gene Ontology (GO) terms associated to each of the obtained hits and returns an evaluated GO annotation for all query sequence(s). Enzyme codes are obtained by mapping to equivalent GOs and InterPro motifs can directly be added to the BLAST based annotation. A basic annotation process with Blast2GO consists of 4 steps (G¨otz et al., 2008):

1. Blasting
2. Interpro-scan
3. GO mapping
4. Functional annotation 

 The Blast2GO Command Line needs a properties file, that contains all the information of the different paramaters that can be changed for the analysis. The properties file can be created with this command:

```bash
blast2go_cli.run -createproperties cli.prop
```

Then use the current version of necessary databases dowloaded, unpackage and imported in the cluster.

> /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/blast2go/DB  

 Replace the cli.prop file in the follow lines with the cluster path described above

 ```
Dbacces.assocdbdata=http://archive.geneontology.org/latest-full/go_201402-assocdb-data.gz
Dbacces.geneinfo=ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz
Dbacces.gene2accession=ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2accession.gz
Dbacces.idmapping=ftp://ftp.pir.georgetown.edu/databases/idmapping/idmapping.tb.gz

 ```

Before start use blast2go_cli, lets process the uniprot aligment:

```bash
DB = /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Trinotate_db/
cp $DB/uniprot_sprot.pep .
 
```

Remember format the database parsing the sequence IDs, because they are needed for the mapping step in Blast2GO:

> The sequences have to be in fasta format and the accession IDs in between “|”. Use `sed -i 's/\s/|/g' uniprot_sprot.pep` to fix the separator line

 

```bash
makeblastdb -dbtype prot -in uniprot_sprot.pep -parse_seqids -out uniprot_sprot 
```

Now let's blast the query sequences against the formatted database, either using the blastx or blastp in command line.

> From testing query sequences remove extra-headers and subset a sample in order to test this workflow (use [subsample.py](https://github.com/RJEGR/infovis/blob/master/subsample.py) in order to process this step).

```bash
 export PATH=/LUSTRE/apps/bioinformatica/blast2GO/blast2go_cli_v1.3.3:$PATH
```

 

```bash
blastx -db ./DB/uniprot/uniprot_sprot -outfmt 5 -evalue 1e-3 -word_size 3 -show_gis -num_alignments 20 -max_hsps 20 -num_threads 24 -out local_blast.xml -query Trinity.fasta.subset 
```

> As authors said, When blasting your sequences, make sure you use the parameter `-show_gis` in order to retrieve the accessions IDs from the formatted database. The resulting file will be an xml file (-outfmt 5), which can be easily loaded in Blast2GO. 

Finally load your blast xml-file results as documentation in blast2go Command Line. 

> Due to licence requierements within the cluster fist start extra-session in:

 

```bash
salloc -N 1 -n 24 --exclusive --nodelist nodo30
ssh nodo30
# password:  ***

```

After finish your task, please exit the node 30 by typing:

```bash
disown %1
exitless
```

> In case to thread the task, please use `nohup` to manage any bast2go task. Ex `nohup blast2go_cli.run command & ` 

## Convert sequences to protein and save as fasta file.

 ```bash
blast2go_cli.run -properties cli.prop -useobo go_latest.obo -loadfasta ./Trinity.fasta.subset -savelorf ./Trinity.fasta.subset.pep

 ```

The follow code Load a DNA fasta file, add the corresponding blast results and perform mapping and annotation. Furthermore, we want to save the .dat file and the PDF report at the current directory with the given name (example).

```bash
blast2go_cli.run -properties cli.prop -useobo DB/go_latest.obo -loadfasta ./Trinity.fasta.subset -loadblast ./local_blast.xml -mapping -annotation -savedat example -savereport example

```

 ## Load b2g file and run Go slim analysis

```bash 
blast2go_cli.run -properties cli2.prop -useobo ../DB/go_latest.obo -loadb2g example_data/example.b2g -goslim example_data/goslim_plant.obo -saveb2g -tempfolder ./tmp
```

## Load fasta file, add corresponding blast result and execute GO mapping and annotation. 

>  runing in nohup ... still running

```bash
nohup \
blast2go_cli.run -properties cli.prop -loadfasta \
Trinity.fasta.subset -loadblast31 \
local_blast.xml -mapping -annotation \
-workspace work_dir -nameprefix prefix -saveb2g -saveannot -savereport \
-saveseqtable -statistics gdatadispie,aecdis -tempfolder ./tmp \
&

```

## Reporte

El programa blast2go_cli.run ubicado dentro del cluster en la ruta PATH parece funcionar con signos de advertencia en los módulos de java (Ej. *java.io.IOException: No space left on device |or| Java HotSpot(TM) 64-Bit Server VM warning: Insufficient space for shared memory file*), ignoró la modalidad en la que el programa  blast2go solicite recursos del cluster, pero reconozco que este demanda memoria para almacenar información temporal, sin embargo el nodo30 (en el cual está instalada la licencia de este programa) está limitado en la memoria para almacenar archivos temporales. 

Por ello, blast2go permite seleccionar otra ruta para  almacenar información temporal (opcion -tempfolder de blast2go_cli.run). Ojo: si no se establece este parámetro, la información de entrada (archivos fasta)  no es cargada satisfactoriamente, y un error de tipo There are no sequences with sequence information se produce sin importar cual sea el modulo del programa blas2go_cli.run que se este tratando de implementar. 

Por default en el manual no utilizan esta configuración temporal y es algo que ya Silvia había anticipado que utilizaramos, para encarar los errores con java sin embargo he notado que es necesario utilizar en todo momento.Por tanto, me di la libertad de probar la herramienta en cuestión con los datos de ejemplo que vienen en la carpeta del PATH/example_data (un subset de datos pequeños de prueba) y FUNCIONA de manera apropiada la mayoría de los módulos que nos interesan implementar (ejemplo 1,5,6,  de la sección 4 del manual); sin embargo, el ejemplo 8 del manual es uno que nos interesa mucho ejecutar por lo que ejecute una prueba el dia Viernes que hasta el dia de hoy NO concluye o modifica el archivo nohup.log desde entonces por lo que decidí abortar la tarea. Cabe mencionar que estas pruebas las hice también implementando datos 'reales' que procesamos en el laboratorio (archivo Trinity.fasta.subset) y todo parece funcionar bien a excepción de la prueba del ejemplo 8 que hasta el momento está en ejecución en el nodo30 en el comando de pid 17007:

```java
java -Xmx100568m -XX:+UseConcMarkSweepGC -XX:+Disable...
```

En la sintaxis: 

```bash

blast2go_cli.run -properties cli.prop -loadfasta \Trinity.fasta.subset -loadblast \local_blast.xml -mapping -annotation -workspace work_dir -nameprefix prefix -saveb2g -saveannot -savereport -saveseqtable -statistics gdatadispie,aecdis -tempfolder ./tmp 
```

Me intriga mucho porque demora tanto la corrida del programa (el mismo comando de java era el que se ejecutaba desde el Viernes con el set de datos de ejemplo), y sospecho por el comando del pid 17007 que el problema del retraso del programa es una vez más con la colector de memoria basura de java. En la idealidad, la herramienta en cuestión disminuiría los tiempos de análisis computacionales que se suelen hacer con su parte gráfica blast2go (visualizar resultados y elaborar algunas inferencias estadísticas), en el caso del archivo Trinity.fasta.subset, es un set de datos muy pequeños (100 secuencias de ADN) que debería ser rápido de analizar.

Sylvia, espero que los esfuerzos de este reporte responda algunas preguntas que en correos anteriores me solicitabas, en caso de que puedas implementar los ejemplos como se describe en el manual los datos de ejemplo se localizan en  PATH/example_data habilitando la configuracion temporal  -tempfolder.Debido a que blast2go_cli es una herramienta de paga, es imperativo que le echemos a correr antes de considerar el comunicarnos con el soporte técnico del software para asistencia (support@blast2go.com).

Quedo pendiente de cualquier pregunta o sugerencia.



>  PATH = /LUSTRE/apps/bioinformatica/blast2GO/blast2go_cli_v1.3.3/
>
> > **Manual**: https://www.blast2go.com/images/b2g_pdfs/blast2go_cli_manual.pdf 
>
> WORKDIR=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/blast2go

## second report

After install blast2go_cli (activated licence) in a cluster node with Linux Red-Hat Enterprise 6.7 system (storage of 131 TB (lustre + infiniband)) + 24 cores) We processed the cli_prop file as manual said: `blast2go_cli.run -createproperties cli.prop` and replace the database path file (downloaded manually) in section  `DataAccessParameters` and start trying the examples described in the manual.Our intention were run GO Mapping,  Annotation  and report blast2go stats at beginning. Therefore, We included the parsing the sequence IDs within the makeblast step plus the fixed reference database (changing separator line with `sed -i 's/\s/|/g' uniprot_sprot.pep`). Then, we ensured to use the parameter `-show_gis` in blastx (2.4.0+) in order to retrieve the accessions IDs from the formatted database. Unfortunately we had unknown of problems trying to test the any example (as the manual described) with either,  the data_example data sets and subset of our own datasets.

The chunk-code of the blastx was as above:

`makeblastdb -dbtype prot -in uniprot_sprot.pep -parse_seqids -out uniprot_sprot`

`blastx -db ./DB/uniprot/uniprot_sprot -outfmt 5 -evalue 1e-3 -word_size 3 -show_gis -num_alignments 20 -max_hsps 20 -num_threads 24 -out local_blast.xml -query Trinity.fasta.subset`

 In the next section I summarize the errors found: *Error "there are not sequences with mapping*" in every use-case-example command. Specially the example 1,2,3 and 5 from section 4.1 in the manual. It issue were tested after the blog indications [here](https://www.blast2go.com/support/blog/22-blast2goblog/111-format-fasta-file-blast) and [Frequently Asked Questions](https://www.blast2go.com/support/faq#faqnoanchor) -No mapping results after loading my own blast xml file without any solution. 

A quick solution to skip this error message were add the -tempfolder option in all the blast2go_cli command nevertheless, only converting sequences to proteins use-case-example finish properly and the rest of the Blast2go_cli use-case-examples stucks running indefinitely without any warning or error message using either, the data_example data sets and subset of our own datasets. Finally we could solve this unreported error by switching the option `-loadblast` to `-loadblast31`. 

The final problem were running the code below than output statistics in pdf file but unfortunately with plots are empty and the log file recall an annotation error: There are no sequences with Mapping, please do the mapping before run Annotation.

`blast2go_cli.run -properties cli.prop -loadfasta example_data/1000_plant.fasta -loadblast31 example_data/1000_plant_blastResult.xml -mapping -annotation -saveb2g example.b2g tempfolder ./tmp -savereport example.pdf`

I attached to this mail the report example.pdf, cli.prop and the log file produced by this step in addition to the own subset of datasets implemented during testing blas2go_cli.

We are very concerned about the time we have had spent using this tool without satisfactory results so we decided to emigrate and use the comprehensive annotation suite trinotate in our lab.