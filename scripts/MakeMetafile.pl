#!/usr/bin/perl

use strict;

open OUT, ">/data/results/$ARGV[0]/$ARGV[0].metafile.txt";

print OUT $ARGV[0], "\n";
print OUT "/data/results/\n";
print OUT "/data/database/GRGENEECnet_2014.txt\n";
print OUT "0.958108\n";
print OUT "0.288343\n";
print OUT "0.782456\n";
print OUT "0\n";
print OUT "0.808274\n";
print OUT "0.886837\n";
print OUT "-1.75286\n";
print OUT "0.1179\n";

close OUT;