#!/usr/bin/perl

use strict;
use warnings;


open(OUT, ">to_run.sh");
print OUT "#!/bin/bash\n";

open(IN, "/gscuser/mwyczalk/projects/CPTAC3/import.CPTAC3b1/BamMap/CPTAC3.b1.RNA-Seq.BamMap.dat");
while(<IN>)
{
	chomp;
	next if $_=~m/^#/;
	my @l = split(/\t/,);
	my $sample=$l[2]."\_\_".substr($l[0], 0, -13);
	system("mkdir -p $sample");
	system("ln -s $l[5] $sample/$sample\_1.fastq.gz");
	my $nxt = <IN>;
	my @m = split(/\t/,$nxt);
        system("ln -s $m[5] $sample/$sample\_2.fastq.gz");
	print OUT "lsf_submit 50 4 $sample.log bash /gscmnt/gc2521/dinglab/qgao/Scripts/RNA/rna_pipeline.sh $sample $sample\_1.fastq.gz $sample\_2.fastq.gz\n";
}

