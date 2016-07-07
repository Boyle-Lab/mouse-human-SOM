#!/usr/bin/perl

package Tree;
require Exporter;

use strict;
use warnings;
use GenericFuncs;

our @ISA = ("Exporter");
our @EXPORT = ("parsenode", "print_nodes");
our $VERSION = 1.00;

sub new {
    return bless {}, 'Tree';
}

sub print_nodes {
    my ($self) = @_;
    my ($left, $right, $dist) = @$self;
    if (ref($left) eq "ARRAY" && ref($right) eq "ARRAY") {
	print "@$right\n";
	print "@$left\n";
        print_nodes($left);
        print_nodes($right);
    } elsif ($left eq "leaf") {
        # we are a leaf in the format ['node', 'node name', distance]
	print "@$self\n";
    }
}

sub infer_ancestral {
    # Infer the ancestral character at the root node by parsimony.
    my ($self, $char_idx) = @_;

#    print STDERR "@$self\n";

    my @chars = ();

    my ($left, $right, $dist) = @$self;
    if ($left eq "leaf") {
	# leaf node
	push @chars, ${$char_idx}{$right};
	return(@chars);

    } elsif (ref($left) eq "ARRAY" && ref($right) eq "ARRAY") {
	
	my @child_chars;
	push @child_chars, infer_ancestral($left, $char_idx);
	push @child_chars, infer_ancestral($right, $char_idx);
#	print STDERR "@$self\t@child_chars\n";

	# Compare the current children with the content of the master
	# chars list. If there is sufficient data to resolve the ancestral
	# base, replace the chars list with the most likely base. Otherwise
	# set the chars list to the union of all children.	
	my $c = get_max_char(\@child_chars);
	if ($c eq "N") {
	    @chars = GenericFuncs::array_union(\@chars, \@child_chars);
	} else {
	    @chars = ($c);
	}
#	print STDERR "@$self\t@chars\n";
    }

    if (!defined($$self[3])) {
	# Root node -- choose the character with highest frequency as the
	# most likely ancestral character and return it (returns "N" if
	# not resolvable, given the data).
	return get_max_char(\@chars);
    } else {
	return @chars;
    }
}

sub get_max_char {
    # Count occurences of characters in an array and return the character
    # with the highest count, or "N" if two or more characters have equal
    # counts.
    my ($chars) = @_;

    my %count;
    foreach my $c (@$chars) {
	$count{$c}++;
    }
    my $max_char;
    my $max = 0;
    foreach my $c (keys(%count)) {
	if ($count{$c} > $max) {
	    $max = $count{$c};
	    $max_char = $c;
	} elsif ($count{$c} == $max) {
	    $max_char = "N";
	}
    }

    return $max_char;
}

sub uproot {
    # Drop the root and right child from a tree, rerooting the tree at the
    # first internal node above.
    my ($self) = @_;
    
    my ($left, $right, $dist, $parent) = @$self;
#    print STDERR "$$left[0]\t$$right[0]\n";
    $$left[3] = undef; # Undefine the parent node
    
    return bless $left, 'Tree';
}
