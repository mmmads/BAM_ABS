#!/usr/bin/perl
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
#This program is to convert ambiguous alignment from SAM format to BED format
#Input: each multiread and its multiple alignments (in SAM format)
#Output: BED format
#BED format: Col1: Chrom, Col2: start posititon, Col3: end, Col4: multiread ids, Col5: flag
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
use warnings;
use strict;
use Data::Dumper;
use Storable;

#Read in file line by line
my $file = "$ARGV[0]"; #Input file
unless(open(IN, $file)) {
    print "Could not open file $file!\n";
    exit;
}
print "Processing $file...\n";
open(OUT,">$file.bed"); #output

while(my $line = <IN>)
{
	my @array = split(/\t/,$line);
	#print $array[10];
	my $end = $array[3] + length($array[9]) - 1; #Compute end position
	print OUT "$array[2]\t$array[3]\t$end\t$array[0]\t$array[1]\n"; 
}
