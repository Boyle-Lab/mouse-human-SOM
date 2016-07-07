#!/usr/bin/perl

package ParseNewick;
require Exporter;
require Tree;

use strict;
use warnings;

our @ISA = ("Exporter");
our @EXPORT = ();
our $VERSION = 1.00;


sub new {
    return bless {}, 'ParseNewick';
}

sub FromStringWithDistance {
    my ($self, $in, $parent)= @_;

    my @node;

    # do we have a distance. If so we want to trim that off and save it for later
    my $distance="";
    my $rparen=rindex($in, ")");
    my $rcolon=rindex($in, ":");
    if ($rparen < $rcolon) {
	$distance = substr($in, $rcolon+1);
	$in = substr($in, 0, $rcolon); # adjust the string to exclude the distance
    }

    # if we don't have a comma we are just the node
    if (index($in, ",") == -1) {
	return ["leaf", $in, $distance, $parent];
    }

    # when we get here, we need to find the first level comma to split on
    my $p = 0;
    # find the split point where we are at one level above nothing!
    if (index($in, "(") != -1) {
	my $depth = 0; # the depth we are in the tree
	for (my $i=0; $i < length($in); $i++) {
	    my $c = substr($in, $i, 1);
	    ++$depth if ($c eq '(');
	    --$depth if ($c eq ')');
	    if ($c eq ',' && $depth == 1) {
		$p = $i;
		last;
	    }
	}
    } else {
	$p=index($in, ",");
    }

    # now we have the split points we need to figure out if the left and right halves begin and end with parens, and ignore them if they do
    my ($start, $endlength) = (0, length($in)-$p-1);
    if (index($in, "(") == 0 && rindex($in, ")") == length($in)-1) {
	$start=1; $endlength--;
    }

    my $left_child = $self->FromStringWithDistance(substr($in, $start, $p-1), \@node);
    my $right_child = $self->FromStringWithDistance(substr($in, $p+1, $endlength), \@node);
    
    @node = ($left_child, $right_child, $distance, $parent);

#    print STDERR "$in\t@node\n";

    return \@node;
}

sub parse {
    my ($self, $tree_str, $node_idx) = @_;
    # remove the trailing ;
    $tree_str =~ s/;s*$//;
    my $root = $self->FromStringWithDistance($tree_str, $node_idx);
    return bless $root, 'Tree';
}
