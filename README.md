# CPTAC3-RNA-related-pipeline

####################
###### Transcript 
####################

# Reference

Genome sequence and annotation are downloaded from Ensembl (GRCh37.75)

# Tools

MapSplice is downloaded from http://www.netlab.uky.edu/p/bioinfo/MapSplice2
Cufflinks is downloaded from http://cole-trapnell-lab.github.io/cufflinks/install/

# Processing

first run `perl get_link.pl /gscuser/mwyczalk/projects/CPTAC3/import.CPTAC3b1/BamMap/CPTAC3.b1.RNA-Seq.BamMap.dat` to generate 'to_run.sh' fil
then run 'bash to_run.sh', which 1) maps fastq file to human genome using MapSplice; 2) processes bam file from MapSplice to estimate transcript expression using Cufflinks; 3) convert GTF file from Cufflinks to BED12 using custom Perl script Convert_GTF_To_Bed12.pl

# Output

Output files can be found in 'TRANSCRIPT_BED/' of each sample folder, reported in standard BED12 format with gene expression incorported in the 4th column (TranscriptID|FPKM)
Detailed about BED12 format can be found in https://genome.ucsc.edu/FAQ/FAQformat.html.

####################
###### Fusion
####################

#Three tools were used for fusion calling
#Their references (GRCh19) were also downloaded from the corresponding database

STAR-Fusion is downloaded from https://github.com/STAR-Fusion/STAR-Fusion/wiki

EricScript is downloaded from https://sites.google.com/site/bioericscript

Integrate is downloaded from https://sourceforge.net/p/integrate-fusion/wiki/Home/

# Fusion calling

first run 'perl get_link.pl /gscuser/mwyczalk/projects/CPTAC3/import.CPTAC3b1/BamMap/CPTAC3.b1.RNA-Seq.BamMap.dat' to generate 'to_run.sh' file

then run 'bash to_run.sh', which runs these three tools individually for each sample

finally run 'perl combine_call.pl DIR (/gscmnt/gc2521/dinglab/qgao/RNA/Batch_20171110)' to merge all the raw fusion calls into one file 'Total_Fusions.tsv'


# Fusion filtering

Since raw fusion calls contain many false positives, extensive filtering was performed.

The basic idea for filtering is:

first get the fusions 1) reported by at least 2 callers, or 2) reported by STAR-Fusion (showing higher sensitivity) but with higher supporting evidence (defined by fusion fragments per million total reads, or FFPM, >0.1)

then remove the fusions in the filtering database located in '/gscmnt/gc2521/dinglab/qgao/Reference/GRCh37.75/FusionDatabase/FilterDatabase', including:
1) uncharacterized genes, immunoglobin genes, mitochondrial genes, etc.
2) fusions from the same gene or paralogue genes (downloaded from https://www.genenames.org/cgi-bin/statistics)
3) fusions reported in normal samples - TCGA normals (from pancan fusion analysis, under review), GTEx tissues (reported in star-fusion output), non-cancer cell study (PMID: 26837576)

In general, simply run 'perl filter.pl' to generate the filtered fusion calls 'Filtered_Fusions.tsv'

Optionally, run 'bash generate_fusion_per_sample.sh' to split the filtered fusion calls into one file per sample.

# Output
In the output file, each row represents one fusion.
There are 9 columns for each fusion:
1) FusionName
2) LeftBreakpoint
3) RightBreakpoint
4) Cancer__Sample
5) JunctionReadCount
6) SpanningFragCount
7) FFPM                 - fusion fragments per million total reads, 'NA' means the fusion was found by both EricScript and Integrate but not STAR-Fusion
8) PROT_FUSION_TYPE     - INFRAME, FRAMESHIFT or '.'
9) CallerN              - number of callers
