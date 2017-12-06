#!/bin/bash

sample=$1
fq_1=$2
fq_2=$3
cpu=4

cd ${sample}

##Examine quality metrics for sequencing reads using fastqc
#mkdir -p FASTQC
#fastqc -o FASTQC -f fastq -t 1 $fq_1 $fq_2

##Mapping
#TopHat
#bowtie2_ref=/gscmnt/gc2521/dinglab/qgao/Reference/GRCh37.75/Bowtie2/hg19
#gtf=/gscmnt/gc2521/dinglab/qgao/Reference/GRCh37.75/GTF/Homo_sapiens.GRCh37.75.gtf
#mkdir -p TOPHAT
#tophat -o TOPHAT -p $cpu -G $gtf $bowtie2_ref $fq_1 $fq_2

#MapSplice
#mkdir -p MAPSPLICE
#bowtie_ref=/gscmnt/gc2521/dinglab/qgao/Reference/GRCh37.75/Bowtie/hg19
#chr_dir=/gscmnt/gc2521/dinglab/qgao/Reference/GRCh37.75/Chromosome
#gtf=/gscmnt/gc2521/dinglab/qgao/Reference/GRCh37.75/GTF/Homo_sapiens.GRCh37.75.gtf
#zcat $fq_1 > MAPSPLICE/left.fastq
#zcat $fq_2 > MAPSPLICE/right.fastq
#python /gscmnt/gc2521/dinglab/qgao/Tools/MapSplice-v2.2.1/mapsplice.py -p $cpu -o MAPSPLICE --bam --gene-gtf $gtf -c $chr_dir -x $bowtie_ref -1 MAPSPLICE/left.fastq -2 MAPSPLICE/right.fastq
#rm -f MAPSPLICE/left.fastq MAPSPLICE/right.fastq
#samtools sort -m 20G MAPSPLICE/alignments.bam MAPSPLICE/sorted.alignments

##Cufflinks
#TopHat bam
#gtf=/gscmnt/gc2521/dinglab/qgao/Reference/GRCh37.75/GTF/Homo_sapiens.GRCh37.75.gtf
#cufflinks -o CUFFLINKS_wTOPHAT -p $cpu --library-type fr-firststrand -g $gtf TOPHAT/accepted_hits.bam

#MapSplice bam
#gtf=/gscmnt/gc2521/dinglab/qgao/Reference/GRCh37.75/GTF/Homo_sapiens.GRCh37.75.gtf
#cufflinks -o CUFFLINKS_wMAPSPLICE -p $cpu --library-type fr-firststrand -g $gtf MAPSPLICE/sorted.alignments.bam

##Generate BED for expressed transcripts
#binDir=/gscmnt/gc2521/dinglab/qgao/Scripts/RNA
#mkdir -p TRANSCRIPT_BED
#perl $binDir/Convert_GTF_To_Bed12.pl CUFFLINKS_wMAPSPLICE/transcripts.gtf exon TRANSCRIPT_BED/${sample}.bed

##Fusion calling
#STAR-Fusion
#genome_lib_dir=/gscmnt/gc2521/dinglab/qgao/Reference/GRCh37.75/FusionDatabase/GRCh37_gencode_v19_CTAT_lib_July192017/ctat_genome_lib_build_dir
#mkdir -p STAR_FUSION
#STAR-Fusion --left_fq $fq_1 --right_fq $fq_2 --CPU $cpu --annotate --examine_coding_effect --genome_lib_dir $genome_lib_dir --output_dir STAR_FUSION

#EricScript
#genome_db=/gscmnt/gc2521/dinglab/qgao/Reference/GRCh37.75/FusionDatabase/ericscript_db_homosapiens_ensembl73
#rm -fr ERICSCRIPT
#ericscript.pl -o ERICSCRIPT --remove -ntrim 0 --refid homo_sapiens -db $genome_db -p $cpu -name $sample $fq_1 $fq_2

#Integrate
#bwts=/gscmnt/gc2521/dinglab/qgao/Reference/GRCh37.75/FusionDatabase/Integrate/bwts
#bam_dir=/gscmnt/gc2521/dinglab/qgao/RNA/Batch_20171110/$sample/TOPHAT/
#fasta=/gscmnt/gc2521/dinglab/qgao/Reference/GRCh37.75/FusionDatabase/Integrate/hg19.fa
#annot=/gscmnt/gc2521/dinglab/qgao/Reference/GRCh37.75/FusionDatabase/Integrate/annot.ensembl.GRCh37.txt
#rm -fr Integrate
#mkdir -p Integrate
#samtools index $bam_dir/accepted_hits.bam
#samtools index $bam_dir/unmapped.bam
#Integrate fusion -reads Integrate/reads.txt -sum Integrate/summary.tsv -ex Integrate/exons.tsv -bk Integrate/breakpoints.tsv -vcf Integrate/bk_sv.vcf -bedpe Integrate/fusions.bedpe $fasta $annot $bwts $bam_dir/accepted_hits.bam $bam_dir/unmapped.bam
