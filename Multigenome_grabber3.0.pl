#!/usr/bin/perl

#This program pulls out promoters of specified length upstream of most of the gene 
#in the 80 genomes of Arabidopsis
#Version 3.0 uses a hash of arrays constructed from the annotation file in order to 
#pull out genomic features as they are being read in from the 80 multigenome file

#Author: David Lowry

use warnings;
use strict;

my @matrix_rows;
my @matrix_array;
my $counter1=0;
my $counter2=0;
my $counter3=0;
my $switch=0;
my @annotation;
my @genes;
my $z=0;
my %Pos_hash;

#Argument 1 is annotation file
my $file1 = $ARGV[0];
open(MATRIX, $file1) or die;

#Argument 3 is the length upstream of the promoter to be extracted
my $length = $ARGV[1];

#Argument 3 is the multigenome dataset
my $file2 = $ARGV[2];
open(GENOMES, $file2) or die;

#Argument 4 is the chromosome, because full genome can be too large memory of hash of arrays
my $chrome = $ARGV[3];

while (<MATRIX>){
	chomp;
	@matrix_rows = split('\t');
	push(@matrix_array, [@matrix_rows]);
	$counter1++;
	}
close MATRIX;

#Create a matrix with all the information for position of promoters
for (my $y=0; $y < $counter1; ++$y){
	if("$matrix_array[$y][2]" eq "gene" && "$matrix_array[$y][6]" eq "+" && "$matrix_array[$y][0]" == "$chrome"){
			my $end=("$matrix_array[$y][3]"-1);
			my $start=($end - $length);
			my $chromosome="$matrix_array[$y][0]";
			my $gene = "$matrix_array[$y][9]";
			my @annotation_rows = ();	
			push (@annotation_rows, $chromosome);
			push (@annotation_rows, $start);
			push (@annotation_rows, $gene);
			push (@annotation, [@annotation_rows]);
			$counter2++;
		}elsif("$matrix_array[$y][2]" eq "gene" && "$matrix_array[$y][6]" eq "-" && "$matrix_array[$y][0]" == "$chrome"){
			my $start="$matrix_array[$y][4]";
			my $chromosome="$matrix_array[$y][0]";
			my $gene = "$matrix_array[$y][9]";
			my @annotation_rows = ();	
			push (@annotation_rows, $chromosome);
			push (@annotation_rows, $start);
			push (@annotation_rows, $gene);
			push (@annotation, [@annotation_rows]);
			$counter2++;
 		}else{
 			next;
 		}
 	}

#Populate hash of arrays with "chromosome\tposition" as the key and gene names as values in arrays
for (my $x=0; $x < $counter2; ++$x){
	my $chr = "$annotation[$x][0]";
	my $start = "$annotation[$x][1]";
	my $gene = "$annotation[$x][2]";
	for (my $z=1; $z < $length + 1; ++$z){
		my $it = $start + $z;
		my $key = "$chr\t$it";
		if ($Pos_hash{$key}[0]){
			push @{ $Pos_hash{$key} }, "$gene";
		}else{
			$Pos_hash{$key} = [$gene];
		}
	}
}

#Read in genomes file and print out promoters by each base pair
while (<GENOMES>){
	chomp $_;
	my ($chr, $pos, @values) = split("\t", $_);
	my $length1 = scalar @values;
	my $find= "$chr\t$pos";
	if ($Pos_hash{$find}[0]){
		my $length2 = scalar @{$Pos_hash{$find}};
		for (my $x=0; $x < $length2; ++$x){
			print "$Pos_hash{$find}[$x]\t$find\t";
			for (my $y=0; $y < $length1-1; ++$y){
				print "$values[$y]\t";
			}
			print "$values[$length1-1]\n";
		}
	}else{
		next;
	}
}