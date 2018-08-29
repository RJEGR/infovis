

# Single rsem sample analysis

Both, the de novo transcriptome assembly and all pre-processed libraries is going to be used as input to perform a sample-specific expression analysis. All reads need to be aligned back against the indexed de novo transcriptome assembled using Bowtie2 (Langmead and Salzberg, 2012), followed by calculation of gene and isoform expression levels using Expectation-Maximization algorithm embedded in the Trinity differential expression modules within Trinity (`align_and_estimate_abundance.pl`) on a per sample basis.

```bash
$UTILS/align_and_estimate_abundance.pl \
                --transcripts Trinity.fasta --seqType fq  \
				--est_method RSEM \
				--aln_method bowtie2 \
                --prep_reference \
                --trinity_mode 	\
                --samples_file samples.file \
				--thread_count=24
```

> The sample file is modified with the current identifier name (Ex. P.qtrim.gz prefix used within trinity –trimmomatic config
>
> Also use tool paths:
>
> BOWTIE2=/LUSTRE/apps/bioinformatica/bowtie2
> UTILS=$OWN/trinityrnaseq/util
> RSEM=/LUSTRE/bioinformatica_data/genomica_funcional/bin/RSEM
> SAM2LS=$OWN/samtools-1.3.1

Final step is output the isoform/genes results in a matrix than could be input in follow statement:

```bash
ls `path_where_rsem_results_are`/rsem/*/*isoforms.results > isoforms.results
```

Then, convert into a matrix as follow and start further analysis and filtering:

```bash
abundance_estimates_to_matrix.pl \
        --gene_trans_map Trinity.gene_trans_map \
        --est_method RSEM \
        --out_prefix RSEM.isoforms \
        --quant_files isoforms.results.txt \
        --name_sample_by_basedir
```

> The `abundance_estimates_to_matrix.pl` script may be included in the trinity/utils folder (Ex. $INSTALLATION_PATH/trinityrnaseq-Trinity-v2.5.1/util/)

> The gene_trans_map file details the relation between gene and isoforms in the assembly within two columns `<tab>` delimited file. It can be done by run get_Trinity_gene_to_trans_map.pl script in `INSTALLATION_PATH/trinity_utils/support_scripts/get_Trinity_gene_to_trans_map.pl trinity.fasta > Trinity.gene_trans_map`



## Visualization

> Ref: https://github.com/bli25broad/RSEM_tutorial

After we ran rsem we can visualize some good plots than model statistics RSEM learned from the data. The resulting file, `namefile.pdf` contains plots of learned fragment length distribution, read length distribution, read start position distribution, quality score information and alignment statistics. For example, the below figure shows the read start position distribution (RSPD) learned from a paired-end data.

Lets suppose we are in the previous work directory (where you did run the `align_and_estimate_abundance.pl` script) and we were imput a sample of name sample1

```bash
cd sample1
RSEM=/LUSTRE/bioinformatica_data/genomica_funcional/bin/RSEM
$RSEM/rsem-plot-model RSEM sample1.stat.pdf

# or
# after abundance step, we can figure-out some figures 
cd T0_R

mv RSEM.stat/ T0_R.stat
cd T0_R.stat

for i in *; do mv $i ${i/RSEM/T0_R}; done

# for better compressive treatment download by copy:
for i in $(find -type d -links 3); do cp -r $i ./download/; done

cd download
for i in $(ls -r); do file=./$i/RSEM.stat/; mv $file ${file/RSEM/$i} ; done

# then run hand-by-hand the script (example):
rsem-plot-model RSEM T3_S4.stat.pdf
```

The fist argument needed in the rsem-plot-model is the prefix name of the sample name.



Let's figure out we're in the rsem folder after

The quality score plot (shown below) plots the observed quality against the theoretical quality score for each nucleotide. In general, sequencing error decreases as the quality score increases. However, the theoretical quality score is more optimistic than the observed quality learned from the data.



figure here



The alignment statistics plot (shown below) is another interesting figure to look at. In this plot, unalignable reads, uniquely mapped reads, multi-mapping reads and reads filtered due to too many alignments are represented by colors of green, blue, grey and red. The x-axis of the histogram bins reads by their number of transcript alignments. The y-axis counts the number of reads belonging to each bin. We also have a pie chart at the top-right corner showing the percentage of each type of read in the data. Based on the pie chart, there are much more multi-mapping reads (X%) than uniquely mapped reads (Y%).



figure here