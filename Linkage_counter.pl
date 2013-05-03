#!/usr/bin/perl

#Author: David Lowry

#This script takes the output of eQTL_linkage_combiner.pl
#and uses it to count the number of local and distant linkages

use warnings;
use strict;


my $counter1=0;
my @link_rows;
my @link_array;
my $counter2=0;
my $counter3=0;
my $counter5=0;



#Open the input file
open LINK, "</Users/dbl2/Documents/Arabidopsis_eQTL/cis_trans_May_5_9_12/Linkages/Combined_pvalue_0.05_May_wet_SNP.csv" or die;


while (<LINK>){
	chomp;
	@link_rows = split(',');
	push(@link_array, [@link_rows]);
	$counter1++;
	}
	
close LINK;

#print "$link_array[2][6]\n";

my $threshold = 0.013;

for (my $x=1; $x < $counter1; ++$x){
	if ($link_array[$x][6] < $threshold){
			$counter5++;
	}else{
		#print "$link_array[$x][0]\n";
	}
}

for (my $x=1; $x < $counter1; ++$x){
	if ($link_array[$x][6] < $threshold && $link_array[$x][1] == $link_array[$x][8]){
		$counter2++;
		if ($link_array[$x][10] > $link_array[$x][3] && $link_array[$x][10] < $link_array[$x][4]){
			#print "$link_array[$x][0]\,$link_array[$x][8]\,$link_array[$x][10]\,local\n";
			$counter3++;
		}else{
			#print "$link_array[$x][0]\,$link_array[$x][8]\,$link_array[$x][10]\,distant\n";
			next;	
		}
	}elsif($link_array[$x][6] < $threshold && $link_array[$x][1]){
		#print "$link_array[$x][0]\,$link_array[$x][8]\,$link_array[$x][10]\,distant\n";
	}else{
		next;
	}
}	

my $all = ($counter1-1); 

my $distant=($counter5-$counter3);

my $per_distant=(($distant/$counter5)*100);

my $per_close=(($counter3/$counter5)*100);

print "Total number of eQTLs of eQTLs examined: $all\n";

print "Total number of eQTLs at P < $threshold: $counter5\n";

print "Total number of eQTLs on the same chromosome: $counter2\n";
	
print "Total number of close linkages: $counter3 ($per_close%)\n";

print "Total number of distant linkages: $distant ($per_distant%)\n";