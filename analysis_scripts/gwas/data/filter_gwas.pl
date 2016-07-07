#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use GenericFuncs;

my $usage = "
filter_gwas.pl

Given two beds, filter bed 1 to remove lines in bed 2 (comparison by name col)

USAGE

compare_gwas.pl <bed 1> <bed 2>

<bed 1>
    Base file.

<query bed>
    File with lines to remove from file 1.

OPTIONS

--help
    Show this message.
\n";

my $ref_file = $ARGV[0];
my $query_file = $ARGV[1];

my $help = 0;

GetOptions (
    "help" => \$help
    );

if ($help || $#ARGV+1 < 2) {
    die "$usage";
}

my %query_features = GenericFuncs::read_into_hash($query_file, 3);

open my $INFILE, '<', $ref_file || die "Cannot read $query_file: $!\n";

while (<$INFILE>) {
    chomp;
    my @tmp = split /\t/, $_;

    if (!exists($query_features{$tmp[3]})) {
	print "$_\n";
    }
}

close $INFILE;
