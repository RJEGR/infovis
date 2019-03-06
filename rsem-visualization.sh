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

rsem-plot-model T3_S4 T3_S4.stat.pdf

# also rename  ############
cd T0_R
mv bowtie2.bam T0_R.transcript.bam

# for i in $(find -type d -links 3); mv $i/bowtie2.bam $i/{$i}.transcript.bam; done


srun rsem-plot-transcript-wiggles --gene-list --show-unique T0_S1 ../Trinity.fasta.gene_trans_map T0_S1_wiggle.pdf &


mv bowtie2.bam T2_R.transcript.bam

# Defining path from trimommatic
OWN=/LUSTRE/bioinformatica_data/RNA/ricardo/bioinformatics
TRIMMOMATIC=$OWN/Trimmomatic-0.36
TRUSEQ=~/
#Defining paths for RSEM
BOWTIE2=/LUSTRE/apps/bioinformatica/bowtie2
UTILS=$OWN/trinityrnaseq/util
RSEM=/LUSTRE/bioinformatica_data/genomica_funcional/bin/RSEM
SAM2LS=$OWN/samtools-1.3.1

#exporting paths, JICase
export PATH=$UTILS:$PATH
export PATH=$BOWTIE2:$PATH
export PATH=$RSEM:$PATH
export PATH=$SAM2LS:$PATH