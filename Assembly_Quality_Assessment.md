# Transcriptome Assembly Quality Assessment

The aim of (denovo) transcriptome assembly is to accurately reconstruct the complete set of transcripts that are represented in the read data (in the absence of reference genome)

In contrast of genome assembly, (**context**: if we have an organism with C chromosomes our optimal (ideally) genome assembly would consist of C long contigs) transcrimptome assembly (optimal) wil vary consist of all posible transcript from all expressed genes. ie. alternatively spliced variants (isoforms). 

Nevertheless there are several contributing factors than negatively affect the accuracy of the construction assembly process.

These factors include error in sequencing process (usually the batch effect is found here), incomplete coverage of transcripts (due to insuficient sequencing depth), real biological variability (ex. alternatively spliced variants) and algorithmic simplification. 

Usualy the most highly expressed transcript do not neccessary constitute the longjest one and the majority of transcripts in a transcriptome assembly will normally have relatively low expression levels.



# N50

Imagine that you line up all the contigs in your assembly in the order of their sequence lengths (Fig. 1a). You have the longest contig first, then the second longest, and so on with the shortest ones in the end. Then you start adding up the lengths of all contigs from the beginning, so you take the longest contig + the second longest + the third longest and so on — all the way until you’ve reached the number that is making up 50% of your total assembly length. That length of the contig that you stopped counting at, this will be your N50 number (Fig. 1b).

We strongly advise against using regular N50 metrics for transcriptome assemblies. Instead, other more appropriate measures can be used. The developers of the transcriptome assembler Trinity have invented [the **ExN50** metric](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats), which takes into account the expression levels of each contig and is therefore a more suitable contig length metric for transcriptomes.  

# Ex90

...

# Reads represented

To assess the read composition of our assembly, we want to capture and count all reads that map to our assembled transcripts, including the properly paired and those that are not we run the process below. Bowtie2 is used to align the reads to the transcriptome and then we count the number of proper pairs and improper or orphan read alignments. First, build a bowtie2 index for the transcriptome and then perform the aligment to catpure the paired-reads aligment statistic. (Ref [here](https://github.com/trinityrnaseq/trinityrnaseq/wiki/RNA-Seq-Read-Representation-by-Trinity-Assembly))



```bash


bowtie2-build Trinity.fasta Trinity.fasta
Then perform the alignment (example for paired-end reads) to just capture the read alignment statistics.

srun bowtie2 -p 10 -q --no-unal -k 20 -x Trinity.fasta -1 R1.P.qtrim.fq -2 R2.P.qtrim.fq | samtools view -@10 -Sb -o ./bowtie2.bam


```

Finally lets summary the bam file using bamtools:

```bash
bamtools stats -in reads_represented/bowtie2.bam
```

> in case you neet, fist install bamtools as follow: https://github.com/pezmaster31/bamtools/wiki/Building-and-installing

A output in your screen will be printed as follow:

**********************************************
Stats for BAM file(s):
**********************************************

Total reads:       13970909
Mapped reads:      13970909	(100%)
Forward strand:    6885283	(49.283%)
Reverse strand:    7085626	(50.717%)
Failed QC:         0	(0%)
Duplicates:        0	(0%)
Paired-end reads:  13970909	(100%)
'Proper-pairs':    11803698	(84.4877%)
Both pairs mapped: 13127153	(93.9606%)
Read 1:            7037955
Read 2:            6932954
Singletons:        843756	(6.03938%)

**********************************************

The [Integrative Genomics Viewer](http://software.broadinstitute.org/software/igv/) is useful for visualizing read support across any of the Trinity assemblies. The bowtie2 alignments generated above, which are currently sorted by read name, can be re-sorted according to coordinate, indexed, and then viewed along with the Trinity assemblies using the IGV browser as follows.

```bash
# sort the alignments by coordinate
samtools sort bowtie2.bam -o bowtie2.coordSorted.bam

# index the coordinate-sorted bam file
samtools index bowtie2.coordSorted.bam

# index the Trinity.fasta file
samtools faidx Trinity.fasta

# view the aligned reads along the Trinity assembly reference contigs.
# note, you can do this by using the various graphical menu options in IGV (load genome 'Trinity.fasta', load file 'bowtie2.coordSorted.bam'), or you can use the command-line tool like so:

igv.sh -g `pwd`/Trinity.fasta  `pwd`/bowtie2.coordSorted.bam
```

