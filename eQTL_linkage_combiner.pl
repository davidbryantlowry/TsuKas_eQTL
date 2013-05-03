#!/usr/bin/perl

#Author: David Lowry

#This script combines the parsed output from low_pvalue_finder2.0.pl
#and combines it with the estimated cM location of genes in the genome
#in order to conduct cis/trans analyses

use warnings;
use strict;

my @gene_rows;
my @gene_array;
my $counter1=0;
my @qtl_rows;
my @qtl_array;
my $counter2=0;



#Open the four files
open GENE, "</Users/dbl2/Documents/Arabidopsis_eQTL/cis_trans_Tom_new_McKay_map_2_1_12/gene_locations_new_McKay_map.csv" or die;
open QTL, "</Users/dbl2/Documents/Arabidopsis_eQTL/cis_trans_May_5_9_12/eQTL_P0.05_wet_May_SNP.csv" or die;


while (<GENE>){
	chomp;
	@gene_rows = split(',');
	push(@gene_array, [@gene_rows]);
	$counter1++;
	}

while (<QTL>){
	chomp;
	@qtl_rows = split(',');
	push(@qtl_array, [@qtl_rows]);
	$counter2++;
	}
	
close GENE;
close QTL;

#print "$qtl_array[1][0]\n";
#print "$gene_array[1][0]\n";
print "gene\,chromosome\,position\,ci_low\,ci_high\,lod\,pvalue\,Gene\,Chr\,Loc\,cM\n";

for (my $x=0; $x < $counter2; ++$x){
	for (my $y=0; $y < $counter1; ++$y){
		if ("$qtl_array[$x][0]" eq "$gene_array[$y][0]"){
				print "$qtl_array[$x][0]\,$qtl_array[$x][1]\,$qtl_array[$x][2]\,$qtl_array[$x][3]\,$qtl_array[$x][4]\,$qtl_array[$x][5]\,$qtl_array[$x][6]\,$gene_array[$y][0]\,$gene_array[$y][3]\,$gene_array[$y][4]\,$gene_array[$y][5]\,$gene_array[$y][6]\n";
			}else{
				next;				
			}
		}
	}