#!/usr/bin/perl

#Author: David Lowry

#This script parses a R/qtl output file in the TabByCol format
#It returns the lowest P-value for each phenotype for FDR analysis 
#in Q-value (Dabney and Storey)

use warnings;
use strict;

my @pval_rows;
my @pval_array;
my $counter1;

my $file = $ARGV[0];
open(PVAL, $file) or die;

#Read lines of R/qtl file into an array of arrays
while (<PVAL>){
	@pval_rows = split(' ');
	$counter1 = scalar @pval_rows;
	push(@pval_array, [@pval_rows]);
	}


for (my $x=1; $x < (($counter1/6)+1); ++$x){
	print"$pval_array[0][($x*6)-1]\,";
	my @ranker = ("$pval_array[1][($x*6)-1]", "$pval_array[2][($x*6)-1]", "$pval_array[3][($x*6)-1]", "$pval_array[4][($x*6)-1]", "$pval_array[5][($x*6)-1]");
	@ranker = sort @ranker;
	#  print "$pval_array[0][($x*6)-1]\t";
	print "$ranker[0]\n";
}

