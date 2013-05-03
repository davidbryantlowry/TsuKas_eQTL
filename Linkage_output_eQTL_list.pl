#!/usr/bin/perl

#Author: David Lowry

#This script takes the output of eQTL_linkage_combiner.pl
#and uses it to output local and distant linkage eQTLs

use warnings;
use strict;


my $counter1=0;
my @link_rows;
my @link_array;
my $counter2=0;
my $counter3=0;
my $counter5=0;



#Open the four files
open LINK, "</Users/dbl2/Documents/Arabidopsis_eQTL/cis_trans_May_5_9_12/Linkages/Combined_pvalue_0.05_May_wet_SNP_hotspot.csv" or die;


while (<LINK>){
	chomp;
	@link_rows = split(',');
	push(@link_array, [@link_rows]);
	$counter1++;
	}
	
close LINK;

#print "$link_array[2][6]\n";

my $threshold = 0.004;

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
			print "$link_array[$x][0]\,$link_array[$x][1]\,$link_array[$x][2]\,$link_array[$x][3]\,$link_array[$x][4]\,$link_array[$x][5]\,$link_array[$x][6]\,$link_array[$x][7]\,$link_array[$x][8]\,$link_array[$x][9]\,$link_array[$x][10]\,$link_array[$x][11]\,local\n";
			$counter3++;
		}else{
			print "$link_array[$x][0]\,$link_array[$x][1]\,$link_array[$x][2]\,$link_array[$x][3]\,$link_array[$x][4]\,$link_array[$x][5]\,$link_array[$x][6]\,$link_array[$x][7]\,$link_array[$x][8]\,$link_array[$x][9]\,$link_array[$x][10]\,$link_array[$x][11]\,distant\n";
			next;	
		}
	}elsif($link_array[$x][6] < $threshold && $link_array[$x][1]){
		print "$link_array[$x][0]\,$link_array[$x][1]\,$link_array[$x][2]\,$link_array[$x][3]\,$link_array[$x][4]\,$link_array[$x][5]\,$link_array[$x][6]\,$link_array[$x][7]\,$link_array[$x][8]\,$link_array[$x][9]\,$link_array[$x][10]\,$link_array[$x][11]\,distant\n";
	}else{
		next;
	}
}	

