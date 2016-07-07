#!/usr/bin/perl

# Given a set of peaks annotated with presence/absence in human, mouse and
# a set of outgroup species, predict whether elements are gained/lost in
# human/mouse based on a simple ancestral reconstruction.

use strict;
use warnings;
use Getopt::Long;
use ParseNewick;
use Tree;
use GenericFuncs;

my $usage = "
infer_gain-loss.pl

For a set of _unmapped_ peaks for human and mouse, compile data on whether
each peak is present or absent in a set of outgroups (assumed to be horse,
dog and elephant), and use the information to infer whether the element
was gained or lost, relative to the ancestral sequence, and on which branch
the event occured.

Peaks are labeled by appending an integer classifier: 0 \= unresolvable,
\(1,2 \= reserved for orthologous peaks\), 3 \= mouse gain, 4 \= human gain,
5 \= mouse loss, 6 = human loss. Ouput sent to STDOUT.

Assumes that all reference peaks are UNmappable between human and mouse!

process_outgroups.pl <unmapped.txt> <mapped_to_horse.txt>
    <mapped_to_dog.txt> <mapped_to_elephant.txt> [OPTIONS]

All input files are in a modified bed format, with the first column
containing the species identifier, followed by the standard bed fields.

<unmapped.txt>
    Reference peaks: all UNMAPPED peaks from human and mouse. Peaks in this
    file labeled as mouse are assumed to have no orthologous sequence in
    the human genome and vice-versa.

<mapped_to_horse.txt> <mapped_to_dog.txt> <mapped_to_elephant.txt>
    Files containing mappings of the peaks in <unmapped.txt> to the
    given genome.

OPTIONS:

--help
    Show this message.    
\n";

my $help = 0;
GetOptions (
    "help" => \$help
    );

if ($help || $#ARGV+1 < 4) {
    die "$usage\n";
}

# Only need mapped features for each outgroup.
my ($ref_f, $equ_f, $can_f, $lox_f, $tree_str) = @ARGV;

# Parse the Newick tree
my $parser = new ParseNewick;
my $tree = $parser->parse($tree_str);

# Prune off a subtree directly above the root for resolving ambiguity
my $pruned = $parser->parse($tree_str);
$pruned = $pruned->uproot();

# Key input hashes on the peak name
my %ref_feats = GenericFuncs::read_into_hash($ref_f, 4);
my %equ_feats = GenericFuncs::read_into_hash($equ_f, 4);
my %can_feats = GenericFuncs::read_into_hash($can_f, 4);
my %lox_feats = GenericFuncs::read_into_hash($lox_f, 4);


foreach my $key (keys(%ref_feats)) {
    my @out = @{$ref_feats{$key}};
    
    my %chars;
    if (${$ref_feats{$key}}[0] eq "hg19") {
	# present in human but not mouse
        $chars{hg19} = 1;
	$chars{mm9} = 0;
    } else {
	# present in mouse but not human
	$chars{hg19} = 0;
        $chars{mm9} = 1;
    }

    if (exists($equ_feats{$key})) {
	$chars{equCab2} = 1;
    } else {
	$chars{equCab2} = 0;
    }
    if (exists($can_feats{$key})) {
	$chars{canFam2} = 1;
    } else {
	$chars{canFam2} = 0;
    }
    if (exists($lox_feats{$key})) {
	$chars{loxAfr3} = 1;
    } else {
        $chars{loxAfr3} = 0;
    }

    my $anc_char = $tree->infer_ancestral(\%chars);

    if ($anc_char eq "N") {
	# Unresolvable ancestral state This means that loxAfr3 is saying
	# something different than equCab2 and canFam2. These can be
	# resolved by dropping loxAfr3 from the tree.
	$anc_char = $pruned->infer_ancestral(\%chars);
    }

    if ($anc_char eq "N") {
	# This really shouldn't happen, but to be safe...
	push @out, 0;
    } elsif ($anc_char == 1) {
	# Element is present at root -- losses
	if (${$ref_feats{$key}}[0] eq "hg19") {
	    # Mouse loss
	    push @out, 5;
	} else {
	    # Human loss
	    push @out, 6;
	}
    } elsif ($anc_char == 0) {
	# Element is not present at root -- gains
	if (${$ref_feats{$key}}[0] eq "hg19") {
            # Human gain
            push @out, 4;
        } else {
            # Mouse gain
            push @out, 3;
        }
    } else { # This should never happen!
	die "Bad character returned by infer_ancestral. Check that all data is valid!\n";
    }

    GenericFuncs::print_array(\@out);

}
