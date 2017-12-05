#!/usr/bin/perl

use strict;
use warnings;

die "Usage: $0  \"gtf\"  \"region\"  \"output basename\" " if (@ARGV < 3);

open (DATA, "$ARGV[0]");
my (%block_count, %strand, %starts, %ends, %chr);
while(<DATA>)
{
	chomp;
	next if $_=~m/#/;
	my @line = split(/\t/,);
	if ($line[2] eq "$ARGV[1]")
	{
		my @ids = split /\"/, $line[8];
		my $trans_id = $ids[3]."\|".$ids[7];
		$block_count{$trans_id}++;
		my $tmp_start = $line[3]-1;
		my $tmp_end = $line[4];
		$strand{$trans_id}=$line[6];
		$starts{$trans_id}.="$tmp_start ";
		$ends{$trans_id}.="$tmp_end ";
		$chr{$trans_id} = "$line[0]";
	}
}

open (OUT, ">$ARGV[2]");
foreach my $id (keys %block_count)
{
	my @lens = ();
	my @block_start = ();
	my @start = split(/ /, $starts{$id});
	my @end = split(/ /, $ends{$id});
	my @sorted_start = sort {$a <=> $b} @start;
	my @sorted_end = sort {$a <=> $b} @end;
	for (my $i=0; $i<=$#start; $i++)
	{
		my $len = $sorted_end[$i]-$sorted_start[$i];
		push @lens, $len;
		my $re_start = $sorted_start[$i]-$sorted_start[0];
		push @block_start, $re_start;
	}
	my $size = join ",", @lens;
	my $block_starts = join ",", @block_start;
	print OUT "$chr{$id}\t$sorted_start[0]\t$sorted_end[-1]\t$id\t0\t$strand{$id}\t$sorted_start[0]\t$sorted_start[0]\t0\t$block_count{$id}\t$size\t$block_starts\n";
}
