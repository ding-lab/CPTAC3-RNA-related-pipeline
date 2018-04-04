#!/usr/bin/perl

use strict;
use warnings;

die "Usage perl filter.pl SAMPLE \n" if (@ARGV < 1);

my $database_dir="/gscmnt/gc2521/dinglab/qgao/Reference/GRCh37.75/FusionDatabase/FilterDatabase";

my (%black, %paralog, %manual, %noncancer, %tcga_normal, %known, %gtex, %large, %tmp);

open(BLACK, "$database_dir/blacklist");
while(<BLACK>)
{
	chomp;
	$black{$_}="";
}

open(PARALOG, "$database_dir/paralog_clusters.dat");
while(<PARALOG>)
{
	chomp;
	my @l=split(/\t/,);
	next if scalar @l < 2;
	for(my $i=0;$i<=$#l;$i++)
	{
		for(my $j=0;$j<=$#l;$j++)
		{
			next if $i==$j;
			my $pp=join("\-\-",$l[$i], $l[$j]);
			$paralog{$pp}="";
		}
	}
}

open(DAT, "$database_dir/blacklist.manual");
while(<DAT>)
{
        chomp;
        $manual{$_}="";
}

open(DAT, "$database_dir/noncancer_cell.txt");
while(<DAT>)
{
        chomp;
        $noncancer{$_}="";
}

open(DAT, "$database_dir/tcga.normal.id.txt");
while(<DAT>)
{
        chomp;
        $tcga_normal{$_}="";
}

#fusions should stay
open(DA, "$database_dir/tcga.published");
while(<DA>)
{
	chomp;
	next if $_ eq "BMPR1B--PDLIM5";
	next if $_ eq "ZC3H7A--BCAR4";
	$known{$_}="";
}

my $dir = $ARGV[0];
my $sample = $ARGV[1];
open(DATA, "$dir/Total_Fusions_in_$sample.tsv");
<DATA>;
open(OUT, ">$dir/Filtered_Fusions_in_$sample.tsv");
while(<DATA>)
{
	chomp;
	my @line = split(/\t/,);
	next if $line[6] eq "NA";
	if(!exists $large{$line[0]})
	{
		$large{$line[0]} = $line[6];
	}else
	{
		if($large{$line[0]} < $line[6])
		{
			$large{$line[0]} = $line[6];
		}
	}
}

open(DATA, "$dir/Total_Fusions_in_$sample.tsv");
my $header = <DATA>;
my @h = split(/\t/, $header);
print OUT join("\t", @h[0..7], $h[9]);

while(<DATA>)
{
        chomp;
        my @line = split(/\t/,);
	my @genes = split(/\-\-/,$line[0]);
	next if $genes[0] eq $genes[1];
#1.fusion from same gene
	next if (exists $black{$genes[0]} || exists $black{$genes[1]});
#2.blacklist
	next if exists $paralog{$line[0]};
#3.fusion from paralog genes
	next if exists $noncancer{$line[0]};
#4.non-cancer fusion
	next if exists $gtex{$line[0]};
#5.gtex fusion
	next if exists $tcga_normal{$line[0]};
#6.tcga normal
	next if exists $manual{$line[0]};
#7.blacklist pairs
	if($line[6] eq "NA")
	{
		print OUT join("\t", @line[0..7], $line[9])."\n";
	}else
	{
		if($large{$line[0]} == $line[6])
		{
			if(!exists $tmp{$line[0]})
			{
				print OUT join("\t", @line[0..7], $line[9])."\n";
			}
			$tmp{$line[0]} = "";
		}
	}
#8. only report breakpoint with the highest score
}
