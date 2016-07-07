#!/usr/bin/perl

package GenericFuncs;
require Exporter;

use strict;
use warnings;

our @ISA = ("Exporter");
our @EXPORT = ("print_array_tab_delim", "print_array_comma_delim",
	       "print_array_space_delim", "read_into_hash");
our $VERSION = 1.00;

my $TMP_F = (time ^ $$ ^ unpack "%32L*", 'ps wwaxl | gzip') . '.tmp';

sub print_array {
    # Print an array as a (tab)-delimited string
    my ($array, $delim, $fh) = @_;

    if (!defined($delim)) {
	$delim = "\t";
    }

    if (defined($fh)) {
        select $fh;
    } else {
        select STDOUT;
    }
    
    my $i;
    for ($i = 0; $i < $#{$array}; $i++) {
        print "${$array}[$i]$delim";
    }
    print "${$array}[$#{$array}]\n";

    if (defined($fh)) {
        select STDOUT;
    }
    return 0;
}

sub print_hash_of_arrays {
    # Print a hash of arrays to given filehandle with given delimiting
    my ($hash, $delim, $fh) = @_;

    if (!defined($delim)) {
        $delim = "\t";
    }
    if (!defined($fh)) {
        $fh = \*STDOUT;
    }

    if (exists(${$hash}{header})) {
	&print_array(${$hash}{header}, $delim, $fh);
    }

    foreach my $key (sort(keys(%{$hash}))) {
	if ($key eq "header") {
            next;
        }
        &print_array(${$hash}{$key}, $delim, $fh);
    }
}

sub read_into_hash {
    # Read a (tab)-delimited file into a hash of arrays
    my ($file, $key_idx, $delim) = @_;

    if (!defined($delim)) {
	$delim = "\t";
    }

    my %hash;

    open my $INFILE, '<', $file || die "Cannot read $file: $!\n";
    while (<$INFILE>) {
	if ($_ =~ /^#/ || $_ !~ m/\S+/g) {
	    # Skip blank and commented lines
	    next;
        }
        chomp;
        my @tmp = split /$delim/, $_;
        my $key = $tmp[$key_idx];
        $hash{$key} = \@tmp;
    }
    close $INFILE;

    return %hash;
}

sub read_into_hash_multikey {
    # Read a (tab)-delimited file into a hash of arrays with keys
    # constructed as the concatenation of a set of columns as
    # specified in @{$key_idx}
    my ($file, $key_idx, $delim) = @_;
    # $key_idx is a reference to an array holding field indeces from
    # which to build the hash keys

    if (!defined($delim)) {
	$delim = "\t";
    }

    my %hash;

    open my $INFILE, '<', $file || die "Cannot read $file: $!\n";
    while (<$INFILE>) {
	if (/^#/) {
            # Skip commented lines
            next;
	}
        chomp;
	my @tmp = split /\t/;
	
	my @keyvals;
	my $i;
	for ($i = 0; $i <= $#{$key_idx}; $i++) {
	    push @keyvals, $tmp[${$key_idx}[$i]];
	}
	my $key = join '_', @keyvals;

	# Check for duplicates in dataset
#	if (exists($hash{$key})) {
#	    print STDERR "$key: $tmp[4]\n";
#	    print STDERR "${$hash{$key}}[4]\n\n";
#	}
	$hash{$key} = \@tmp;
    }
    close $INFILE;

    return %hash;
}

sub orth_table {
    # Build a two-way orthologs table (as a hash) from a file.
    my ($orth_f) = @_;

    my %orthologs;

    open my $IN, '<', $orth_f;
    while (<$IN>) {
	chomp;
	my @tmp = split /\t/, $_;

	# Expected column format puts human orth in 3, mouse orth in 4
	# and format is species_orth in that column
	my $h_orth = $tmp[3];
	my $m_orth = $tmp[4];
	
	$h_orth =~ s/human_//;
	$m_orth =~ s/mouse_//;

	# Some entries have a suffix included in the gene name: "__N" where
	# N is a number indicating (something??) -- drop this
	$h_orth =~ s/__\d+//;
	$m_orth =~ s/__\d+//;

	$orthologs{$h_orth}{mm} = $m_orth;
	$orthologs{$h_orth}{hg} = $h_orth;
	$orthologs{$m_orth}{mm} = $m_orth;
	$orthologs{$m_orth}{hg} = $h_orth;

    }
    close $IN;

    return %orthologs;
}

sub check_orths {
    # Given a species, (string of comma-delimited) gene name(s) and an
    # orthologs table, get the applicable gene names for human and mouse.
    my ($spp, $keys_str, $orthologs) = @_;
    
    my @keys = split /,/, $keys_str;
    
    my $hg_gene;
    my $mm_gene;
    foreach my $key (@keys) {
        if ($spp eq "hg19") {
            $hg_gene = $key;
#           print STDERR "hg_gene: $key, ";
            if (exists(${$orthologs}{$hg_gene})) {
                $mm_gene = ${$orthologs}{$hg_gene}{mm};
#               print STDERR "mm_gene: $mm_gene\n";
                last;
            }
        } else {
            $mm_gene = $key;
#           print STDERR "mm_gene: $key";
            if (exists(${$orthologs}{$mm_gene})) {
                $hg_gene = ${$orthologs}{$mm_gene}{hg};
#		print STDERR ", hg_gene: $hg_gene\n";
		last;
	    }
	}
    }
    return ($hg_gene, $mm_gene);
}
