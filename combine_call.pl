#!/usr/bin/perl

use strict;
use warnings;

die "Usage perl combine_call.pl DIR \nFilterDB and STARDB are set to be static in the script\n" if (@ARGV < 1);

my $filterDB = "/gscmnt/gc2521/dinglab/qgao/Reference/GRCh37.75/FusionDatabase/FilterDatabase";
open(IN, "$filterDB/hgnc_complete_set.txt");
#different names for the same gene
<IN>;
my %hash;
while(<IN>)
{
        chomp;
        my @l=split(/\t/,);
        $hash{$l[1]}=$l[0];
        if($l[10] ne "")
        {
                $l[10]=~s/\"//g;
                my @m=split(/\|/,$l[10]);
                foreach my $g (@m)
                {
                        $hash{$g}=$l[0];
                }
        }
}

my $dir = $ARGV[0];
my (@names, $sample, $cancer, %mark, %star, %integrate, %eric, %gtex);

foreach my $input (glob("$dir/*__*/STAR_FUSION/star-fusion.fusion_predictions.abridged.annotated.coding_effect.tsv"))
{
	@names = split(/\//,$input);
	$sample = $names[-3];
	open(IN, "$input");
	<IN>;
	while(<IN>)
	{
		chomp;
		my @l = split(/\t/,);
		my $id = join("\_\_", $sample, $l[0]);
		if($_=~m/GTEx_Recurrent/)
		{
			$gtex{$l[0]} = 1;
		}
		$star{$id} = "$l[0]\t$l[5]\t$l[7]\t$sample\t$l[1]\t$l[2]\t$l[13]\t$l[19]";
		my @m = split(/\-\-/, $l[0]);
		if(exists $hash{$m[0]})
		{
			$mark{$hash{$m[0]}}=$m[0];
		}
		if(exists $hash{$m[1]})
		{
			$mark{$hash{$m[1]}}=$m[1];
		}
		#always use the names in STAR-Fusion
	}
}

foreach my $input (glob("$dir/*__*/ERICSCRIPT/*.results.total.tsv"))
{
	@names=split(/\//,$input);
        $sample = $names[-3];
	open(IN, "$input");
	<IN>;
	while(<IN>)
	{
		chomp;
		my @l=split(/\t/,);
		my $a=(exists $hash{$l[0]} && exists $mark{$hash{$l[0]}}) ? $mark{$hash{$l[0]}} : $l[0];
        	my $b=(exists $hash{$l[1]} && exists $mark{$hash{$l[1]}}) ? $mark{$hash{$l[1]}} : $l[1];
		if($l[3]=~m/Unable/)
		{
			$l[3]="Unable";
		}
		if($l[6]=~m/Unable/)
		{
			$l[6]="Unable";
		}
		my $id = join("\_\_", $sample, $a."\-\-".$b);
		$eric{$id} = "$a\-\-$b\tchr$l[2]\:$l[3]\:$l[4]\tchr$l[5]\:$l[6]\:$l[7]\t$sample\t$l[10]\t$l[11]";
#exact breakpoint info will be extracted using Integrate
	}
}

foreach my $input (glob("$dir/*__*/Integrate/summary.tsv"))
{
        @names=split(/\//,$input);
        $sample = $names[-3];
        open(IN, "$input");
        <IN>;
	open(IN2, "$dir/$sample/Integrate/breakpoints.tsv");
	<IN2>;
        while(<IN>)
        {
                chomp;
                my @l=split(/\t/,);
		my $match=<IN2>;
		chomp $match;
		my @ll = split(/\t/,$match);
                my $a=(exists $hash{$l[1]} && exists $mark{$hash{$l[1]}}) ? $mark{$hash{$l[1]}} : $l[1];
                my $b=(exists $hash{$l[2]} && exists $mark{$hash{$l[2]}}) ? $mark{$hash{$l[2]}} : $l[2];
		my $id = join("\_\_", $sample, $a."\-\-".$b);
		$integrate{$id} = "$a\-\-$b\tchr$ll[2]\:$ll[4]\:+\tchr$ll[5]\:$ll[7]\:+\t$sample\t$l[6]\t$l[7]";
#strand info will be corrected using EricScript
        }
}

open(OUT1, ">Total_Fusions.tsv");
open(OUT3, ">STARFusion_not.tsv");
print OUT1 "FusionName\tLeftBreakpoint\tRightBreakpoint\tCancer\_\_Sample\tJunctionReadCount\tSpanningFragCount\tFFPM\tPROT_FUSION_TYPE\tGTEx\tCallerN\n";
print OUT3 "#FusionName\tJunctionReadCount\tSpanningFragCount\tSpliceType\tLeftGene\tLeftBreakpoint\tRightGene\tRightBreakpoint\tLargeAnchorSupport\tLeftBreakDinuc\tLeftBreakEntropy\tRightBreakDinuc\tRightBreakEntropy\tFFPM\n";

my $stardb = "/gscmnt/gc2521/dinglab/qgao/Reference/GRCh37.75/FusionDatabase/GRCh37_gencode_v19_CTAT_lib_July192017/ctat_genome_lib_build_dir";
open(REF, "$stardb/ref_annot.gtf");
my %refdb;
my %transname;
while(<REF>)
{
        chomp;
	my $id;
        next if $_=~m/^#/;
        my @l = split(/\t/,);
        if($l[2] eq "gene")
        {
                my @a=split(/\"/,$l[8]);
                $transname{$a[9]}=$a[1];
#GRC37 and 38 have different formats
        }elsif($l[2] eq "exon")
        {
                $id=$l[0].":".$l[3].":".$l[6];
                $refdb{$id}="";
                $id=$l[0].":".$l[4].":".$l[6];
                $refdb{$id}="";
        }
}

my @tmp_sample;
foreach my $g (keys %eric)
{
	my $beforeformat;
        if(exists $integrate{$g} && !exists $star{$g})
        {
                if($eric{$g}=~m/Unable/)
                {
			my @e = split(/\t/,$eric{$g});
			my @i = split(/\t/,$integrate{$g});
			my ($b1, $b2);
			if($e[1]=~m/Unable/)
			{
				my @a=split(/\:/,$e[1]);
				my @b=split(/\:/,$i[1]);
				$b1=join("\:",$a[0],$b[1],$a[2]);
			}else
			{
				$b1=$e[1];
			}
                        if($e[2]=~m/Unable/)
                        {
                                my @a=split(/\:/,$e[2]);
                                my @b=split(/\:/,$i[2]);
                                $b2=join("\:",$a[0],$b[1],$a[2]);
                        }else
                        {
                                $b2=$e[2];
                        }
			
			$beforeformat = "$e[0]\t$b1\t$b2\t$e[3]\t$e[4]\t$e[5]";
                }else
                {
                        $beforeformat = "$eric{$g}";
                }
		my @before = split(/\t/,$beforeformat);
		my @m = split(/\-\-/,$before[0]);
		if($m[0] eq "PPP2R4")
		{
                	$m[0] = "PTPA";
		}
		#one case with different names
		my $type=(exists $refdb{$before[1]} && exists $refdb{$before[2]}) ? "ONLY_REF_SPLICE" : "INCL_NON_REF_SPLICE";
		print OUT3 "$before[0]\t$before[4]\t$before[5]\t$type\t$m[0]\^$transname{$m[0]}\t$before[1]\t$m[1]\^$transname{$m[1]}\t$before[2]\tYES_LDAS\tGT\t1\tAG\t1\t0.1\n";
		push @tmp_sample, $before[3];
        }
}

my $annodir = "/gscmnt/gc2741/ding/qgao/tools/STAR-Fusion/FusionAnnotator";
my $new_annot = `$annodir/FusionAnnotator --genome_lib_dir $stardb --annotate STARFusion_not.tsv | $annodir/util/fusion_to_coding_region_effect.pl --fusions - --genome_lib_dir $stardb`;
my @lines = split(/\n/, $new_annot);
shift @lines;
for(my $i=0; $i<=$#lines; $i++)
{
        my @after = split(/\t/,$lines[$i]);
	if($lines[$i]=~m/GTEx_Recurrent/)
	{
		$gtex{$after[0]} = 1;
	}
}

foreach my $g (keys %star)
{
	my @l = split(/\t/,$star{$g});
	my $val = (exists $gtex{$l[0]}) ? 1 : 0;
        if(exists $integrate{$g})
        {
                if(exists $eric{$g})
                {
                        print OUT1 "$star{$g}\t$val\t3\n";
                }else
                {
                        print OUT1 "$star{$g}\t$val\t2\n";
                }
        }else
        {
                if(exists $eric{$g})
                {
                        print OUT1 "$star{$g}\t$val\t2\n";
                }else
                {
                        print OUT1 "$star{$g}\t$val\t1\n";
                }
        }
}

for(my $i=0; $i<=$#lines; $i++)
{
	my @after = split(/\t/,$lines[$i]);
	my $val = (exists $gtex{$after[0]}) ? 1 : 0;
	print OUT1 "$after[0]\t$after[5]\t$after[7]\t$tmp_sample[$i]\t$after[1]\t$after[2]\tNA\t$after[19]\t$val\t2\n";
}
system("rm -f STARFusion_not.tsv");
