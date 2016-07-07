#!/usr/bin/perl

# Given a gene expression counts matrix, orthologs list and GC content tables 
# for two species, build a table that associates genes with their respective GC
# content from each species.

use strict;
use warnings;
use GenericFuncs;

my $genes_mtrx_f = $ARGV[0];
my $orth_f = $ARGV[1];
my $annots_1_f = $ARGV[2];
my $annots_2_f = $ARGV[3];

# Read the orthologs table into a pair of hashes
my %orthologs = GenericFuncs::orth_table($orth_f);

# Read the annotations into hashes keyed by gene_id (not quoted)
my %annots_1 = GenericFuncs::read_into_hash($annots_1_f, 0);
my %annots_2 = GenericFuncs::read_into_hash($annots_2_f, 0);

my @header = ("mouse_name", "mouse_GC", "human_name", "human_GC");
GenericFuncs::print_array(\@header);

# Loop over the expression file to produce the GC table. It is important
# to maintain the line-order of the input file!
my $i = 0;
open my $EXPR, '<', $genes_mtrx_f || die "Cannot read $genes_mtrx_f: $!\n";
while (<$EXPR>) {
    # Skip the header line
    if (!$i++) {
	next;
    }

    my @line = split /\t/, $_;

    my ($mm_gene, $hg_gene) = ($line[0], $line[1]);
    
    my ($gc_mm, $gc_hs) = ("NA", "NA");

    if (defined($mm_gene) && $mm_gene ne "NULL" && exists($annots_2{$mm_gene})) {
	$gc_mm = ${$annots_2{$mm_gene}}[2];
    }
    if (defined($hg_gene) && $hg_gene ne "NULL" && exists($annots_1{$hg_gene})) {
        $gc_hs = ${$annots_1{$hg_gene}}[2];
    }

    my @out = ($mm_gene, $gc_mm, $hg_gene, $gc_hs);
    GenericFuncs::print_array(\@out);
}
close $EXPR;
