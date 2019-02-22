#!/usr/bin/perl

use strict;
#use lib "/ifs/scratch/c2b2/dv_lab/plata/NETWORK/";
#use Bio::SeqIO;


#my $seq_in = Bio::SeqIO -> new(-file => $ARGV[0], -format => 'fasta');

#while (my $seq_obj = $seq_in -> next_seq()){
#  my $id = $seq_obj->display_id();
#  print "$id\tchr\t+\t1\t2\t3\t4\n";
#}
my $dir1="/data/genomes/";
my $dir2="/data/results/";
my $genome = "$dir1/$ARGV[0]";
my $outfile = "/data/results/$ARGV[0]/$ARGV[0].dummy_order";

open FASTA, $genome;
open OUT, ">$outfile";
while (<FASTA>){
  if ($_ =~ />(\S+)/){
    my $id=$1;
    print OUT "$id\tchr\t+\t1\t2\t3\t4\n";
  }
}

close FASTA;
close OUT;