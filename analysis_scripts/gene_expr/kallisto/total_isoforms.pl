#!/usr/bin/perl

use strict;
use warnings;
use GenericFuncs;
use Getopt::Long;

my $usage = "

Given a table of gene expression values, possibly containing
multiple rows for isoforms of a given gene, produce a new table
with a single gene value representing the sum of all isoforms.

USAGE:

total_isoforms.pl <infile>

<infile>
    A tab-delimited text file containing expression values in one
    or more columns.

OPTIONS:

--g-col
    Column in <infile> containing identifiers.
    Default = 0

--total-cols
    Columns to sum, as a comma-delimited string. Multiple,
    consecutive columns can be specified as \'n:m\'
    Default = 1:12

--weighted-eff-len
    For each gene+cell, compute abundance-weighted average
    over transcripts as the total effective length.

--help
    Show this message.
\n";

my ($infile) = @ARGV;

my $help = 0;
my $g_col = 0;
my $cols_str = "1:12";
my $weighted_len = 0;
my $len_col = 13;

GetOptions (
    "help" => \$help,
    "g-col=i" => \$g_col,
    "total-cols=s" => \$cols_str,
    "weighted-eff-len" => \$weighted_len
    );

if ($help || $#ARGV+1 < 1) {
    die "$usage\n";
}

# Process the cols string
my @cols;
my @strs = split /,/, $cols_str;
foreach my $str (@strs) {
    my @tmp = split /:/, $str;
    
    my $first = $tmp[0];
    my $last = $tmp[$#tmp];
    for (my $i = $first; $i <= $last; $i++) {
	push @cols, $i;
    }
}

# Go through the file and build the output row by row, totalling rows for isoforms as we find them.
my %dat;
my $rows = 0;
my $dups = 0;

my @header;

open my $IN, '<', $infile || die "Cannot read $infile: $!\n";
while (<$IN>) {
    chomp;
    my @tmp = split /\t/, $_;

    if ($_ =~ m/^#/) {
	push @tmp, "wgt_avg_eff_len";
	@header = @tmp;
	next;
    }

    if (!exists($dat{$tmp[$g_col]})) {
	# New record
	$dat{$tmp[$g_col]}{data}[0] = \@tmp;
	$dat{$tmp[$g_col]}{nrows} = 1;	
	$rows++;
    } else {	
	$dups++;
	push @{$dat{$tmp[$g_col]}{data}}, \@tmp;
	$dat{$tmp[$g_col]}{nrows}++;
    }
}

# Total/average over rows
my %out = average_genes(\%dat, \@cols, $g_col);
$out{header} = \@header;

# Print the result
GenericFuncs::print_hash_of_arrays(\%out, "\t", \*STDOUT);


sub average_genes {
    my ($dat, $cols, $g_col) = @_;

    my %out;    

    foreach my $key (sort(keys(%{$dat}))) {
	my $eff_len = get_avg_eff_len($dat{$key});
	my @row = @{$dat{$key}->{data}->[0]};
	for (my $i = 1; $i < $dat{$key}->{nrows}; $i++) {
	    foreach my $j (@{$cols}) {
#           GenericFuncs::print_array($out{$tmp[$g_col]}, "\t", \*STDERR);                                                                                                                                          
		$row[$j] += $dat{$key}->{data}->[$i]->[$j];
#           GenericFuncs::print_array($out{$tmp[$g_col]}, "\t",\*STDERR);                                                                                                                                           
	    }

	}
	push @row, $eff_len;
	$out{$key} = \@row;
    }

    return %out;
}

sub get_avg_eff_len {
    my ($dat) = @_;

#    print STDERR $dat->{nrows}, "\n";

    my @tmp;
    my $gene_ab = 0;
    for (my $i = 0; $i < $dat->{nrows}; $i++) {
#	GenericFuncs::print_array($dat->{data}->[$i], "\t", \*STDERR);
	my $grp_ab = $dat->{data}->[$i]->[11] + $dat->{data}[$i]->[12];
	$tmp[$i][1] = $grp_ab;
	if ($grp_ab == 0) {	    
	    next;
	}
	$tmp[$i][0] += $dat->{data}->[$i]->[13];
	$gene_ab += $grp_ab;
    }

#    print STDERR "$gene_ab\n";

    if ($gene_ab == 0) {
	return 0;
    }
    my $eff_len = 0;
    for (my $i = 0; $i < $dat->{nrows}; $i++) {
	if ($tmp[$i][1] == 0) {
	    next;
	}
	$eff_len += $tmp[$i][0] * ($tmp[$i][1] / $gene_ab);
    }

    return $eff_len;
}
