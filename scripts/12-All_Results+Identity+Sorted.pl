#!/usr/bin/perl

use strict;

if (scalar @ARGV == 0){
  die "USAGE: $0 <Genome>\n";
}

my $dir1 ="/data/results/$ARGV[0]";
my $dir2 ="/data/database";

my @p=`cat $dir1/$ARGV[0].r1.new.GLOBUS.Gene_Prob.txt`;
my @i=`cat $dir2/GRGENEECnet_2014.txt`;
my @l=`cat $dir1/$ARGV[0].gene`;
my $total = scalar @l * 10000;


my %index;
foreach my $line (@i){
  chomp $line;
  my @parts = split /\t/, $line;
  $index{$parts[0]}=$parts[1];
}
my %iden;
my %candidates;
my %gene_index;

foreach my $line (@l){
  chomp $line;
  my @sections = split /\s+\/\/\/\s+/, $line;
  my @parts = split /\t/, $sections[0];
  my $n = $parts[2];
  my @iparts = split /\t/, $sections[1];
  for (my $count=0; $count<$n; $count++){
    push @{$candidates{$parts[0]}}, $parts[$count+3];
    $gene_index{$parts[0]}=$parts[1];
    my $ec=$index{$parts[$count+3]};
    my $iden=$iparts[$count];
    #print "$parts[1]\t$ec\t$iden\n";
    $iden{$parts[1]}{$ec}=$iden;
  } 
}

my %out;
my %sor;

foreach my $line (@p){
  chomp $line;
  my @parts = split /\t/, $line;
  next if (@parts ==2);
  my $n = $parts[1];
  my %p;
  next unless (exists $candidates{$parts[0]});
  my @candidates = @{$candidates{$parts[0]}};
  for (my $count = 0; $count < $n; $count++){
    $p{$candidates[$count]}=sprintf ("%.3f", $parts[$count+2]/$total);
  }
  my $out = sprintf ("%.3f", $parts[$n+2]/$total);
  my $name = $gene_index{$parts[0]};
  #print "$name";
  foreach my $candidate (sort {$p{$b} <=> $p{$a}} keys %p){
  #foreach my $candidate (sort keys %p){
    $sor{"$name\t$index{$candidate}\t$p{$candidate}"}=$p{$candidate};
    $out{"$name\t$index{$candidate}\t$p{$candidate}"}="$name\t$index{$candidate}\t$p{$candidate}\t$iden{$name}{$index{$candidate}}\n";
    #last; ### Only best P.
  }
  #print "\n";
}


foreach my $pair (sort {$sor{$b} <=> $sor{$a}} keys %sor){
#foreach my $pair (sort keys %sor){
#foreach my $pair (sort keys %sor){
  print $out{$pair};
}




