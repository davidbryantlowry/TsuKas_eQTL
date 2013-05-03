#!/usr/bin/perl

#Author: David Lowry

#Version 2.0 does not return the lowest p-value.  Instead it returns
#info on the location of all QTLs discovered below a particular P-Value
#threshold.  

use warnings;
use strict;

my @pval_rows;
my @pval_array;
my $counter1;


my $file = $ARGV[0];
open(PVAL, $file) or die;

#Read lines of genotype file into an array
while (<PVAL>){
	@pval_rows = split('\s');
	$counter1 = scalar @pval_rows;
	push(@pval_array, [@pval_rows]);
	}
	
print "gene,chromosome,position,ci_low,ci_high,lod,pvalue\n";

for (my $x=1; $x < (($counter1/6)+1); ++$x){
	  for (my $y=1; $y < 6; ++$y){
	  		if (($pval_array[$y][($x*6)-1]) < 0.05){
	  			print "$pval_array[0][($x*6)-1]\,";
	  			print "$pval_array[$y][($x*6)-6]\,";
	  			print "$pval_array[$y][($x*6)-5]\,";
	  			print "$pval_array[$y][($x*6)-4]\,";
	  			print "$pval_array[$y][($x*6)-3]\,";
	  			print "$pval_array[$y][($x*6)-2]\,";
	  			print "$pval_array[$y][($x*6)-1]\n";
	  		}else{
	  			next;
	  	}
	  }
}

