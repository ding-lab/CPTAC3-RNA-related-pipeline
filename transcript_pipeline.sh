#!/bin/bash

sample=$1
fq_1=$2
fq_2=$3
cpu=4

cd ${sample}

#MapSplice
mkdir -p MAPSPLICE
bowtie_ref=/gscmnt/gc2521/dinglab/qgao/Reference/GRCh37.75/Bowtie/hg19
chr_dir=/gscmnt/gc2521/dinglab/qgao/Reference/GRCh37.75/Chromosome
gtf=/gscmnt/gc2521/dinglab/qgao/Reference/GRCh37.75/GTF/Homo_sapiens.GRCh37.75.gtf
zcat $fq_1 > MAPSPLICE/left.fastq
zcat $fq_2 > MAPSPLICE/right.fastq
python /gscmnt/gc2521/dinglab/qgao/Tools/MapSplice-v2.2.1/mapsplice.py -p $cpu -o MAPSPLICE --bam --gene-gtf $gtf -c $chr_dir -x $bowtie_ref -1 MAPSPLICE/left.fastq -2 MAPSPLICE/right.fastq
rm -f MAPSPLICE/left.fastq MAPSPLICE/right.fastq
samtools sort -m 20G MAPSPLICE/alignments.bam MAPSPLICE/sorted.alignments

##Cufflinks
gtf=/gscmnt/gc2521/dinglab/qgao/Reference/GRCh37.75/GTF/Homo_sapiens.GRCh37.75.gtf
cufflinks -o CUFFLINKS_wMAPSPLICE -p $cpu --library-type fr-firststrand -g $gtf MAPSPLICE/sorted.alignments.bam

##Generate BED for expressed transcripts
binDir=/gscmnt/gc2521/dinglab/qgao/Scripts/RNA
mkdir -p TRANSCRIPT_BED
perl $binDir/Convert_GTF_To_Bed12.pl CUFFLINKS_wMAPSPLICE/transcripts.gtf exon TRANSCRIPT_BED/${sample}.bed

