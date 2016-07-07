#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use GenericFuncs;

my $usage = "
remap_gwas.pl

Convert GWAS SNP associtions, in bed format, from one genome assembly to
another by matching the rsid to the dbSNP reference set. Prints mapped
features to STDOUT and unmappable features to STDERR.

To start with GWAS Catalog records directly from EMBL, first convert
into bed format using the following command:

grep -v \"DATE ADDED\" EMBL_gwas.txt | awk -F \$'\t' '{if (\$12 != \"\" && \\
\$13 != \"\") printf \"chr%s\\t%s\\t%s\\t%s\\t%s;%s\\n\", \$12, \$13, \$13+1, \$22, \\
\$35, \$36}' > EMBL_gwas.hg38.bed

Make sure to use the dbSNP build matching the build for the GWAS Catalog
records or you will get a large number of unmapped features. As of
October 2015, this is dbSNP buld 142.

USAGE

remap_gwas.pl <reference bed> <query bed> > mapped.bed 2> unmapped.bed

<reference bed>
    dbSNP catalog, containing all chromosomes.

<query bed>
    EMBL GWAS catalog, in bed format.

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

my %ref_features = GenericFuncs::read_into_hash($ref_file, 3);

open my $INFILE, '<', $query_file || die "Cannot read $query_file: $!\n";

while (<$INFILE>) {
    chomp;
    my @tmp = split /\t/, $_;

    if (exists($ref_features{$tmp[3]})) {
	$tmp[0] = ${$ref_features{$tmp[3]}}[0];
	$tmp[1] = ${$ref_features{$tmp[3]}}[1];
	$tmp[2] = ${$ref_features{$tmp[3]}}[2];
	GenericFuncs::print_array(\@tmp, "\t", \*STDOUT);
    } else {
	print STDERR "$_\n";
    }
}

close $INFILE;
