#!/usr/bin/perl

use strict;

my $dir="/data/results/$ARGV[0]";
my $gene_file = "$dir/$ARGV[0].gene";
my $bbhs = "$dir/$ARGV[0].blast.bbh.out";

open GENE, $gene_file;
my %genes;
while (<GENE>){
  my @parts = split /\t/, $_;
  $genes{$parts[1]}++;
}

close GENE;

open BBH, $bbhs;
open OUT, ">$dir/$ARGV[0].blast.bbh.filtered.out";
while (<BBH>){
  my @parts = split /\t/, $_;
  if (exists $genes{$parts[0]}){
    print OUT $_;
  }
}

close BBH;
close OUT;
