#!/usr/bin/perl

use strict;
use warnings;

die "Usage: $0  \"file_from_Matt\" " if (@ARGV < 1);

open(OUT, ">to_run.sh");
print OUT "#!/bin/bash\n";

open(IN, "$ARGV[0]");
<IN>;
my $header = <IN>;
die "Matt has changed the format of input file!\n" if $header!~m/^#/;
die "Matt has changed the format of input file!\n" if $header!~m/Case/;
die "Matt has changed the format of input file!\n" if $header!~m/Disease/;
die "Matt has changed the format of input file!\n" if $header!~m/DataPath/;
my @h = split(/\t/, $header);
my ($sample_idx, $disease_idx, $datapath_idx);
for(my $i=0; $i<=$#h; $i++)
{
        $sample_idx = $i if $h[$i] eq "Case";
        $disease_idx = $i if $h[$i] eq "Disease";
        $datapath_idx = $i if $h[$i] eq "DataPath";
}

open(IN, "sort -k 1,1 $ARGV[0] |");
while(<IN>)
{
	chomp;
	next if $_=~m/^#/;
	my @l = split(/\t/,);
	my $sample=$l[$disease_idx]."\_\_".$l[$sample_idx];
	my $nxt = <IN>;
	my @m = split(/\t/,$nxt);
	die "Not paired!\n" if $l[$sample_idx] ne $m[$sample_idx];
	system("mkdir -p $sample");
        system("ln -s $l[$datapath_idx] $sample/$sample\_1.fastq.gz");
        system("ln -s $m[$datapath_idx] $sample/$sample\_2.fastq.gz");
	print OUT "lsf_submit 50 4 $sample.log bash /gscmnt/gc2521/dinglab/qgao/Scripts/RNA/rna_pipeline.sh $sample $sample\_1.fastq.gz $sample\_2.fastq.gz\n";
}

