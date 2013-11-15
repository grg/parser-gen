package ParserJSONImport;

use 5.010001;
use strict;
use warnings;

use vars qw(@EXPORT @EXPORT_OK);
use Exporter;

our @ISA = qw(Exporter);
our $VERSION = '0.01';
@EXPORT = qw(readJSONHeaders);
@EXPORT_OK = qw();

use Carp;
use Parser;
use JSON;

# Read a set of headers from a JSON file
sub readJSONHeaders {
	my $filename = shift;

	my $json = JSON->new->allow_nonref;

	my $expandedFilename = $filename;
	$expandedFilename =~ s{ ^ ~ ( [^/]* ) }
		      { $1
			      ? (getpwnam($1))[7]
			      : ( $ENV{HOME} || $ENV{LOGDIR}
				    || (getpwuid($<))[7]
                                )
	  }ex;

	# Open, and read the file
	my $json_text;
	my $size = (stat($expandedFilename))[7] or
		croak "cannot open '$filename' for reading: $!";
	open (my $fh, "<", $expandedFilename) or
		croak "cannot open < $filename: $!";
	read $fh, $json_text, $size;
	close ($fh);

	# Decode the JSON
	my $hdrInfo = $json->decode($json_text);

	# Walk through the structure and extract relevant data
	my $firstHdr = $hdrInfo->{'firstHdr'};
	my $headers = $hdrInfo->{'headers'};

	# Clean up the headers
	foreach my $hdrName (keys(%$headers)) {
		my $hdr = $headers->{$hdrName};

		$hdr->{'enumName'} = enumName($hdrName);

		# Fill in missing fields
		$hdr->{'len'} = 0 if !defined($hdr->{'len'});
		$hdr->{'lenBytes'} = [] if !defined($hdr->{'lenBytes'});
		$hdr->{'lenMap'} = {} if !defined($hdr->{'lenMap'});
		$hdr->{'nxtHdrBytes'} = [] if !defined($hdr->{'nxtHdrBytes'});
		$hdr->{'nxtHdrMap'} = {} if !defined($hdr->{'nxtHdrMap'});

		# Identify valid lengths
		my $lengths = [];
		if ($hdr->{'len'} > 0) {
			$lengths->[0] = $hdr->{'len'};
		} else {
			push @$lengths, sort(values(%{$hdr->{'lenMap'}}));
		}
		$hdr->{'lengths'} = $lengths;

		# Massage the extract fields into something slightly nicer
		my $newExtractFields = [];
		for (my $i = 0; $i < scalar(@{$hdr->{'extractFields'}}); $i += 3) {
			my $fieldName = $hdr->{'extractFields'}->[$i];
			my $fieldPos = $hdr->{'extractFields'}->[$i+1];
			my $fieldLen = $hdr->{'extractFields'}->[$i+2];
			push @$newExtractFields, {
				'name' => $fieldName,
				'pos' => $fieldPos,
				'len' => $fieldLen
			};
		}
		$hdr->{'extractFields'} = $newExtractFields;
	}

	# Clean up the header sequences
	my $hdrSeqs = [];
	for my $hdrSeq (@{$hdrInfo->{'hdrSeqs'}}) {
		my $newHdrSeq = [];
		my %refCounts;
		for (my $i = 0; $i < scalar(@$hdrSeq); $i += 3) {
			my $hdrName = $hdrSeq->[$i];
			my $hdrPos = $hdrSeq->[$i+1];
			my $hdrLen = $hdrSeq->[$i+2];

			my $refCountName = $hdrName;
			my $hdr = $headers->{$hdrName};
			$refCountName = $hdr->{'refCountName'} if defined($hdr->{'refCountName'});
			$refCounts{$refCountName} = 0 if !defined($refCounts{$refCountName});
			my $cnt = ++$refCounts{$refCountName};
			push @$newHdrSeq, {
				'name' => $hdrName,
				'pos' => $hdrPos,
				'len' => $hdrLen,
				'cnt' => $cnt
			};
		}
		push @$hdrSeqs, $newHdrSeq;
	}

	return $headers, $firstHdr, $hdrSeqs;
}


1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

ParserJSONImport - Perl extension for blah blah blah

=head1 SYNOPSIS

  use ParserJSONImport;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for ParserJSONImport, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.


=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

A. U. Thor, E<lt>grg@localdomainE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2012 by A. U. Thor

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.1 or,
at your option, any later version of Perl 5 you may have available.


=cut
