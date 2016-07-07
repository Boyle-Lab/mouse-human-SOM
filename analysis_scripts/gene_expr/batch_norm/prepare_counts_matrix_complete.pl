#!/usr/bin/perl

# Prepare raw counts matrix for normalization in R, starting with featureCounts
# output data. Does not filter based on orthologs table -- places all rows
# with a value in at least one dataset in the output table.

use strict;
use warnings;
use GenericFuncs;

my $infile_GM12878 = $ARGV[0];
my $infile_K562 = $ARGV[1];
my $infile_CH12 = $ARGV[2];
my $infile_MEL = $ARGV[3];
my $orth_f = $ARGV[4];

# Read the counts into hashes keyed by gene name
my %GM12878 = GenericFuncs::read_into_hash($infile_GM12878, 0);
my %K562 = GenericFuncs::read_into_hash($infile_K562, 0);
my %CH12 = GenericFuncs::read_into_hash($infile_CH12, 0);
my %MEL = GenericFuncs::read_into_hash($infile_MEL, 0);

my %orthologs = GenericFuncs::orth_table($orth_f);

my @out;

# Loop over one file for each species, building/updating the given matrix rows\
# as we go
my %genes;
$genes{index} = 0;

foreach my $key (sort(keys(%GM12878))) {
    my @row = @{$GM12878{$key}};

    my ($hg_gene, $mm_gene) =  GenericFuncs::check_orths("hg19", $row[0],
							 \%orthologs);
#    print STDERR "$hg_gene, $mm_gene\n";

    &update_row(\@out, \%genes, \%MEL, \%CH12, \%K562, \%GM12878,
		$hg_gene, $mm_gene);
}


foreach my $key (sort(keys(%CH12))) {
    my @row = @{$CH12{$key}};

    my ($hg_gene, $mm_gene) =  GenericFuncs::check_orths("mm9", $row[0],
							 \%orthologs);
    
#    print STDERR "$hg_gene, $mm_gene\n";
    
    &update_row(\@out, \%genes, \%MEL, \%CH12, \%K562, \%GM12878,
		$hg_gene, $mm_gene);
}

sub update_row {
    my ($out, $genes, $MEL, $CH12, $K562, $GM12878, $hg_gene, $mm_gene)  = @_;

    my $id;
    my $exists = 0;
    if (defined($hg_gene) && exists(${$genes}{$hg_gene})) {
	$id = ${$genes}{$hg_gene};
	$exists = 1;
    } elsif (defined($mm_gene) && exists($genes{$mm_gene})) {
	$id = ${$genes}{$mm_gene};
	$exists = 1;
    } else {
	$id = ${$genes}{index};
	if (defined($hg_gene)) {
	    ${$genes}{$hg_gene} = $id;
	}
	if (defined($mm_gene)) {
	    ${$genes}{$mm_gene} = $id;
	}
	${$genes}{index}++;
    }

#    print STDERR "$id, $mm_gene, $hg_gene\n";

    if (!$exists) {
	# Create a new row at index $id
	my @out_row = ($mm_gene, $hg_gene, 0, 0, 0, 0);
	
	if (defined($mm_gene)) {

	    my (@m_row, @c_row);

	    if (defined(${$MEL}{$mm_gene})) {
		@m_row = @{${$MEL}{$mm_gene}};
	    }
	    if (defined(${$CH12}{$mm_gene})) {
		@c_row = @{${$CH12}{$mm_gene}};
	    }

	    if (!defined($out_row[0])) {
		$out_row[0] = $mm_gene;
	    } 	    
	    
	    if (@m_row) {
		$out_row[2] = $m_row[6] + $m_row[7];
	    }
	    if (@c_row) {
		$out_row[3] = $c_row[6] + $c_row[7];
	    }

	}
	if (defined($hg_gene)) {
	    
	    my (@k_row, @g_row);

	    if (defined(${$K562}{$hg_gene})) {
		@k_row = @{${$K562}{$hg_gene}};
	    }
	    if (defined(${$GM12878}{$hg_gene})) {
		@g_row = @{${$GM12878}{$hg_gene}};
	    }

	    if (!defined($out_row[1])) {
		$out_row[1] = $hg_gene;
            }
	    
	    if (@k_row) {
		$out_row[4] = $k_row[6] + $k_row[7];
	    }
	    if (@g_row) {
		$out_row[5] = $g_row[6] + $g_row[7];
	    }
	}

	$out[$id] = \@out_row;

    } else {
	# Get the existing row and make any updates necessary.
	if (defined($mm_gene)) {

	    if (!defined(${${$out}[$id]}[0])) {
                ${${$out}[$id]}[0] = $mm_gene;
            }
	    
	    my (@m_row, @c_row);

	    if (defined(${$MEL}{$mm_gene})) {
		@m_row = @{${$MEL}{$mm_gene}};
	    }
	    if (defined(${$CH12}{$mm_gene})) {
		@c_row = @{${$CH12}{$mm_gene}};
	    }

	    if (@m_row) {
		${${$out}[$id]}[2] = $m_row[6] + $m_row[7];
	    }
	    if (@c_row) {
		${${$out}[$id]}[3] = $c_row[6] + $c_row[7];
	    }

        }
	if (defined($hg_gene)) {

	    if (!defined(${${$out}[$id]}[1])) {
                ${${$out}[$id]}[1] = $hg_gene;
            }

	    my (@k_row, @g_row);

	    if (defined(${$K562}{$hg_gene})) {
		@k_row = @{${$K562}{$hg_gene}};
	    }
	    if (defined(${$GM12878}{$hg_gene})) {
		@g_row = @{${$GM12878}{$hg_gene}};
	    }

	    if (@k_row) {
		${${$out}[$id]}[4] = $k_row[6] + $k_row[7];
	    }
	    if (@g_row) {
		${${$out}[$id]}[5] = $g_row[6] + $g_row[7];
	    }
	}
    }
}

# Print the result
my @colnames = ("mm_gene", "hg_gene", "MEL", "CH12", "K562", "GM12878");
GenericFuncs::print_array(\@colnames);
foreach my $row (@out) {
    if (!defined($$row[0])) {
	$$row[0] = "NULL";
    }
    if (!defined($$row[1])) {
	$$row[1] = "NULL";
    }
    GenericFuncs::print_array($row);
}
