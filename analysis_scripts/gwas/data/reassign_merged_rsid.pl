#!/usr/bin/perl

# Search ENTREZ using the eutils portal at
# http://eutils.ncbi.nlm.nih.gov/entrez/eutils/
# for dbsnp rsid's and replace the rsid in the input file with the current
# rsid if the given rsid was merged into another one.

use strict;
use warnings;
use WWW::Mechanize;
use Getopt::Long;
use XML::Simple;
use GenericFuncs;

my $usage =
"\n
Searches ENTREZ through the eutils portal
(http://eutils.ncbi.nlm.nih.gov/entrez/eutils/)
to find the active rsid for dbSNP records that have been merged and
update the input file accordingly.

USAGE:

reasssign_merged_rsid.pl <input bed> > out.bed

OPTIONS:

--help
    Show this message.

CREDITS AND LICENSE:

Copyright (C) 2016, Adam Diehl

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Any questions/comments regarding this program may be directed to Adam Diehl:
adadiehl\@umich.edu
\n";

if ($#ARGV + 1 < 1) {
    die "$usage\n";
}

my $infile = $ARGV[0];

my $help = 0;

GetOptions (
    "help" => \$help,
    );

# Check for help option and exit with usage message if found.
if ($help) {
    die "$usage\n";
}

# Set up the virtual browser
my $mech = WWW::Mechanize->new(
    ssl_opts => {
        verify_hostname => 0,
        # Quick and dirty hack to avoid errors about missing SSL certificates.
    },
    );

my $i = 0;
open my $INFILE, '<', $infile || die "Cannot read $infile:$!\n";
while (<$INFILE>) {
    chomp;
    my @tmp = split /\t/, $_;

    my $uid = $tmp[3];
    $uid =~ s/rs//;
    
    # This may need to be updated if it changes in the future!
    my $URL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=snp&id=' . $uid;
    print STDERR "Query URL: $URL\n";
    
    # Run the query and process results
    my %result = &get_result($mech, $URL);

    if (!%result) {
	print STDERR "No result found for: $_\n";
	next;
    }

    if ($uid eq $result{SNP_ID}) {
	print "$_\n";
    } else {
	print STDERR "\tUpdating rsid to rs$result{SNP_ID}\n";
	$tmp[3] = "rs" . $result{SNP_ID};
	GenericFuncs::print_array(\@tmp, "\t", \*STDOUT);
	$i++;
    }
}
close $INFILE;

print STDERR "Done. Updated $i records.\n";

######################
# Subroutines

sub get_result {
    # Retrieve XML data from a URL and return data as a hash.
    my $mech = $_[0];
    my $url = $_[1];

    $mech->get($url);
    
    my $xml = $mech->content;
    my $xml_obj = XMLin($xml);

    my %res;

    # Handle multiple results by returning an empty hash
    if (!(ref($xml_obj->{DocSum}) eq 'HASH')) {
	return %res;
    }

    $res{Id} = $xml_obj->{DocSum}->{Id};

    foreach my $row (@{$xml_obj->{DocSum}->{Item}}) {
	if (exists($row->{content})) {
	    $res{$row->{Name}} = $row->{content};
	}
    }

#    foreach my $key (keys(%res)) {
#	print STDERR "$key\t$res{$key}\n";
#    }

    return %res;
}

sub save_json {
    # Save the JSON data for a file we have downloaded.
    my $json = $_[0];
    my $download_path = $_[1];
    my $out_root = $_[2];

    my $outfile = $download_path . $out_root
	. '.' . ${$json}{accession} . '.json';
    open my $OUT, '>', $outfile;
    print $OUT to_json($json, {pretty=>1});
    close $OUT;
}
