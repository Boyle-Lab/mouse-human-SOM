#!/usr/bin/perl

# Convert gene names from one system to another based on a translation table.

use strict;
use warnings;
use GenericFuncs;
use Getopt::Long;

my $usage = "
Convert gene names in a file from one system to another based on a translation
table. Uses hashes to compare keys from the specified column in the input file to
locate the matching identifier in the tranlation file, replacing the value with
the specified column of the transation table. Results are printed to STDOUT.

Can generically be used to reassign columns in one file based on the contents
of another file.

USAGE:

rename_genes.pl <infile> <trans-file> [OPTIONS]

<infile>
    A tab-delimited text file containing values associated wtih a given set of
    identifiers.

<trans-file>
    A tab-delimited text file containing identifiers matched to those in <infile>,
    along with alternative identifiers (one set per row).

OPTIONS:

--g-col
    Column in <infile> containing identifiers.
    Default = 0

--t-col
    Column in <trans-file> containing matching identifiers.
    Default = 0

--r-col
    Column in <trans-file> from which to take replacement identifiers.
    Default = 1

--help
    Show this message.
\n";

my ($genes_f, $trans_f) = @ARGV;

my $g_col = 0;
my $t_col = 0;
my $r_col = 1;
my $help = 0;

GetOptions (
    "help" => \$help,
    "g-col=i" => \$g_col,
    "r-col=i" => \$r_col,
    "t-col=i" => \$t_col
    );

if ($help || $#ARGV+1 < 2) {
    die "$usage\n";
}

my %genes = GenericFuncs::read_into_hash($genes_f, $g_col);
my %trans = GenericFuncs::read_into_hash($trans_f, $t_col);

my $match = 0;
my $no_match = 0;
foreach my $key (sort(keys(%genes))) {
    my @row = @{$genes{$key}};

    if (exists($trans{$key})) {
	$row[$g_col] = ${$trans{$key}}[$r_col];
	$match++;
#	GenericFuncs::print_array(\@row);
    } else {
#	GenericFuncs::print_array(\@row, "\t", \*STDERR);
	$no_match++;
    }
    GenericFuncs::print_array(\@row);
}

print STDERR "Done: Updated $match rows.\n";

if ($no_match) {
    print STDERR "Warning: Identifiers for $no_match records were not updated because there was no match in $trans_f!\n";
}
