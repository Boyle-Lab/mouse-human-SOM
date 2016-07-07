#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use GenericFuncs;

my $usage = "

\n";

my $genes_f = $ARGV[0];
my $orths_f = $ARGV[1];

my $help = 0;
GetOptions (
    "help" => \$help
    );

if ($help || $#ARGV+1 < 2) {
    die "$usage\n";
}

my @idx = (0, 1);
my %expr = GenericFuncs::read_into_hash_multikey($genes_f, \@idx);
my %orthologs = GenericFuncs::orth_table($orths_f);

my %genes;
$genes{index} = 0;

my @out;

foreach my $key (sort(keys(%expr))) {

    my @row = @{$expr{$key}};

#    print STDERR "$row[0], $row[3]\n";

    my $hg_gene;
    my $mm_gene;
    if ($row[0] eq "hg19") {
	($hg_gene, $mm_gene) = GenericFuncs::check_orths("hg19", $row[1], \%orthologs);
    } else {
	($hg_gene, $mm_gene) = GenericFuncs::check_orths("mm9", $row[1], \%orthologs);
    }

    my $id;
    my $exists = 0;
    if (defined($hg_gene) && exists($genes{$hg_gene})) {
	$id = $genes{$hg_gene};
	$exists = 1;
    } elsif (defined($mm_gene) && exists($genes{$mm_gene})) {
	$id = $genes{$mm_gene};
	$exists = 1;
    } else {
	$id = $genes{index};
	if ($row[0] eq "hg19") {
	    $genes{$hg_gene} = $id;
	} else {
	    $genes{$mm_gene} = $id;
	}
	$genes{index}++;
    }
    
    if (!$exists) {
	# Create a new row at index $id
	my @out_row = ($mm_gene, $hg_gene, 0, 0, 0, 0, 0, 0, 0, 0);

	if ($row[0] eq "mm9") {
	    if (!defined($out_row[0])) {
		$out_row[0] = $mm_gene;
	    }
	    $out_row[2] = $row[13];
	    $out_row[3] = $row[12];
	    $out_row[6] = $row[17];
	    $out_row[7] = $row[16];	    
	} else {
	    if (!defined($out_row[1])) {
		$out_row[1] = $hg_gene;
            }
	    $out_row[4] = $row[13];
            $out_row[5]= $row[12];
            $out_row[8] = $row[17];
            $out_row[9] = $row[16];
	}

	$out[$id] = \@out_row;

    } else {
	# Get the existing row and make any updates necessary.
	if ($row[0] eq "mm9") {
	    if (!defined(${$out[$id]}[0])) {
                ${$out[$id]}[0] = $mm_gene;
            }
	    ${$out[$id]}[2] = $row[13];
            ${$out[$id]}[3]= $row[12];
            ${$out[$id]}[6] = $row[17];
            ${$out[$id]}[7] = $row[16];
        } else {
	    if (!defined(${$out[$id]}[1])) {
                ${$out[$id]}[1] = $hg_gene;
            }
            ${$out[$id]}[4] = $row[13];
            ${$out[$id]}[5]= $row[12];
            ${$out[$id]}[8] = $row[17];
            ${$out[$id]}[9] = $row[16];
	}
    }
}

GenericFuncs::print_array([("#mm_gene", "hg_gene", "fpkm_m", "fpkm_c", "fpkm_k",
			    "fpkm_g", "q_m", "q_c", "q_k", "q_g")]);
foreach my $row (@out) {
    if (!defined($$row[0])) {
	$$row[0] = "NULL";
    }
    if (!defined($$row[1])) {
	$$row[1] = "NULL";
    }
    GenericFuncs::print_array($row);
}

