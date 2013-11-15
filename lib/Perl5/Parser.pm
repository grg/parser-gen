package Parser;

use warnings;
use strict;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);
use Exporter;

@ISA = qw(Exporter);
$VERSION = '1.00';
@EXPORT = qw(setWordWidth adjPos wordStartPos enumName
	calcParserLocs getParserLocs
	getWorstPktInfo getWorstHdrInfo
	getMaxHdrSeqLen getMaxHdrSeq
	maxParserSeqLen
	setHeaders getHeaderNames numHeaders getHeader
	getFirstHdrName getFirstHdrDecPos getFirstHdrEnum
	cmpNum
	sortParserSeqLen
	sortParserSeqAlpha
	genericProcParams
	makeParserName
	makeParserIfaceName
	log2 log10
	genPktByteStream
	genExtractByteStream
	printVerilogArray
	setHdrSeqs getHdrSeqs
	minPktSize
	genHdrPosCheck
	getLastDecByte
	calcOffsetsFromLen printLenCaseBranch
	calcExtractTables getExtractTables
	calcIfcArrays getIfcArrayPos getIfcArrayLen getIfcArrayLenDW
	createIfcInsts
	copyIfcFromArray copyIfcToArray
	getFinalParsers getPrecedingParsers
	getMaxHdrLen getHdrWidth
	fmtWidthVal
	getBufRdWidth fmtBRAmt
	fmtVal
	setShiftWidth getShiftWidth
	calcProgLookupTable getProgLookupTable
	setProgExtractWidth setProgLookupInputs setProgLookupWidth setProgMaxRange setProgStateWidth
	setProgFirstLookupWordsNew getProgFirstLookupWordsNew 
	expandFilename
	parseTCAMEntry parseFirstLookupEntry
	setTCAMTableFile getTCAMTableFile
	);
@EXPORT_OK = qw();

use Carp;
use constant {
	DONE	=> 'DONE',
	FCS_LEN	=> 4,
	IPG	=> 20, # 12 (IPG) + 8 (preamble)
	MIN_PKT_SIZE => 60,
	MAX_HDR_LEN => 255,
	PROG_LOOKUP_DONE => 255
};

# Programmable field names
use constant {
	HIST 		=> 'HIST',
	STATE 		=> 'STATE',
	LOOKUP_STR 	=> 'LOOKUP_STR',
	NEXT_STATE 	=> 'NEXT_STATE',
	LOOKUP_OFFSETS 	=> 'LOOKUP_OFFSETS',
	EXTRACT_OFFSETS => 'EXTRACT_OFFSETS',
	EXTRACT_DSTS	=> 'EXTRACT_DSTS',
	SHIFT 		=> 'SHIFT',
	OFFSET_INDEX 	=> 'OFFSET_INDEX',
	OFFSET_AMT 	=> 'OFFSET_AMT',
	NEXT_LOOKUP_OFFSETS => 'NEXT_LOOKUP_OFFSETS',
};
my @allProgFields = (HIST, STATE, LOOKUP_STR, NEXT_STATE, LOOKUP_OFFSETS,
	EXTRACT_OFFSETS, EXTRACT_DSTS, SHIFT, OFFSET_INDEX, OFFSET_AMT,
	NEXT_LOOKUP_OFFSETS);

use List::Util qw(max min);
use POSIX;


my $wordByteWidth = 0;
my $headers;
my $hdrSeqs;
my $numHeaders = 0;
my %parsers;
my @parserSeqs;
my %wideParsers;
my $maxLag = 0;
my $maxLagSeq;
my $maxHdrSeqLen = 0;
my $maxHdrSeq;
my $worstPktRxCycles = -1;
my $worstPktParseCycles = -1;
my $worstPktParseSeq;
my $worstHdrRxCycles = -1;
my $worstHdrParseCycles = -1;
my $worstHdrParseSeq;

my %extractByteCnt;
my %extractMap;
my %extractHdrCnt;
my $totalExtractBytes = 0;
my $totalExtractEntries = 0;

my $firstHdrName;
my $firstHdrDecPos;
my %parserIfcMap;
my $numIfcSW = 0;
my $numIfcDW = 0;
my $shiftWidth = 0;

my %progTableEntries;
my $progNumEntries = 0;
my $progUniqueEntries = 0;
my $progExtraEntries = 0;
my $progFirstLookupWords;
my $progFirstLookupWordsNew;

my $progExtractWidth = 0;
my $progLookupInputs = 0;
my $progLookupWidth = 0;
my $progMaxRange = 0;
my $progStateWidth = 0;

my $tcamTableFile = '';


# How wide are words?
sub setWordWidth {
	$wordByteWidth = shift;
}

# Calculate the minimum start position to encompas the specified locaiton
sub adjPos {
    my $val = shift;
    return $val >= $wordByteWidth ? $val - $wordByteWidth + 1 : 0;
}

# Calculate the minimum start position to encompas the specified locaiton
sub wordStartPos {
    my $val = shift;
    $val = 0 if $val < 0;
    return $val - ($val % $wordByteWidth);
}

# Convert a header name to an enum name
sub enumName {
	my $hdrName = shift;
	$hdrName =~ s/-/_/g;
	$hdrName =~ s/\+/__/g;
	return $hdrName;
}

# Get the last decision byte for a given header
sub getLastDecByte {
	my $hdrName = shift;
	my $hdr = getHeader($hdrName);

	my $lastDecByte = -1;
	foreach my $byte (@{$hdr->{lenBytes}}) {
		$lastDecByte = $byte if $byte > $lastDecByte;
	}
	foreach my $byte (@{$hdr->{nxtHdrBytes}}) {
		$lastDecByte = $byte if $byte > $lastDecByte;
	}

	return $lastDecByte;
}

# Calculate where parsers should be located
sub calcParserLocs {
	my $hdrSeqs = shift;
	my $allowWide = shift;
	$allowWide = 0 if !defined($allowWide);

	$firstHdrName = undef;
	$firstHdrDecPos = 0;
	%parsers = ();
	@parserSeqs = ();
	%wideParsers = ();
	$maxLag = -1;
	$maxHdrSeqLen = -1;
	$worstPktParseCycles = -1;
	$worstPktRxCycles = -1;
	$worstPktParseSeq = undef;
	my $worstPktParseRatio = 0;
	$worstHdrParseCycles = -1;
	$worstHdrRxCycles = -1;
	$worstHdrParseSeq = undef;
	my $worstHdrParseRatio = 0;
	my %pSeqStrs;
	foreach my $hdrSeq (@$hdrSeqs) {
		my $base = 0;
		my $parserSeq = [];
		my $pSeqStr = "";
		my $lag = 0;
		my $parseCycles = 0;
		my $decLen = 0;
		my $currSeqLen = 0;
		foreach my $hdrInfo (@$hdrSeq) {
			my $hdrName = $hdrInfo->{'name'};
			my $hdrPos = $hdrInfo->{'pos'};
			my $hdrLen = $hdrInfo->{'len'};
			my $hdrLastDecByte = getLastDecByte($hdrName) + 1;

			if ($hdrPos == 0) {
				$firstHdrName = $hdrName;
				$firstHdrDecPos = $hdrLastDecByte;
			}

			my $parseWidth = $wordByteWidth;
			$parseWidth *= 2 if ($allowWide and
				$lag >= $wordByteWidth and $currSeqLen < 2 * $wordByteWidth);

			# Record the current sequence only if:
			#   - there's data in the sequence, AND:
			#     - the sequence occupies exactly parseWidth bytes, OR
			#     - the decision point occupies exactly parseWidth bytes, OR
			#     - the current decision point spills over to the next word
			if ($currSeqLen != 0 and
			    ($currSeqLen % $parseWidth == 0 or
			     $decLen % $parseWidth == 0 or
		             ($currSeqLen > $parseWidth and
		              int($decLen / $wordByteWidth) != int($currSeqLen / $wordByteWidth)) or
		             ($currSeqLen > $parseWidth and
		              int($decLen / $wordByteWidth) == int($currSeqLen / $wordByteWidth) and
			      $currSeqLen % $wordByteWidth + $hdrLastDecByte > $wordByteWidth) or
		             ($currSeqLen < $parseWidth and
			       $currSeqLen + $hdrLastDecByte > $parseWidth))) {

				# Default: increase lag by the unprocessed data in
				# final word with prev header
				# Note: NEVER increases by more than wordByteWidth
				my $lagDelta = (-($hdrPos - $base));
				if ($allowWide and $lag >= $wordByteWidth and
					$currSeqLen <= $parseWidth) {
					$lagDelta += $wordByteWidth;
				} else {
					$lagDelta %= $wordByteWidth;
				}

				# Adjustment of lagDelta if the decision point
				# is *beyond* the end of the header.
				# e.g. MPLS: look at the first byte of nxt hdr
				if ($decLen > $parseWidth and
				    POSIX::ceil($decLen / $wordByteWidth) > POSIX::ceil($currSeqLen / $wordByteWidth)) {
					my $revisit = POSIX::ceil($decLen / $wordByteWidth) * $wordByteWidth;
					$revisit -= POSIX::ceil($currSeqLen / $wordByteWidth) * $wordByteWidth;
					$lagDelta += $revisit;
				}

				# Record the parser sequence
				push (@parserSeqs, $parserSeq) if (!defined($pSeqStrs{$pSeqStr}));
				$pSeqStrs{$pSeqStr} = 1;

				# Update the lag
				# (Different calcs based upon whether we're skipping bytes)
				if ($decLen % $parseWidth == 0 or
				    ($currSeqLen > $parseWidth and
		                     int($decLen / $wordByteWidth) != int($currSeqLen / $wordByteWidth))) {
					my $lastDecWordByte = $decLen + $parseWidth - 1;
					$lastDecWordByte -= $lastDecWordByte % $parseWidth;
					$lagDelta = -($currSeqLen - $lastDecWordByte);
				}
				$lag += $lagDelta;
				$lag %= $wordByteWidth if ($lag < 0);

				# Stats tracking
				if ($lag > $maxLag) {
					$maxLag = $lag;
					$maxLagSeq = $hdrSeq;
				}

				$parserSeq = [];
				$pSeqStr = "";
				$decLen = 0;
				$currSeqLen = 0;
				$base = $hdrPos;

				$parseWidth = $wordByteWidth;
				$parseWidth *= 2 if ($allowWide and $lag > $wordByteWidth);
			}

			my $offset = $hdrPos - $base;
			$offset %= $wordByteWidth if $offset >= $parseWidth;
			$parsers{$offset} = {} if (!defined($parsers{$offset}));
			$parsers{$offset}{$hdrName} = 1;
			if ($allowWide and $lag >= $wordByteWidth and $hdrLastDecByte > $wordByteWidth) {
				$wideParsers{$offset} = {} if (!defined($wideParsers{$offset}));
				$wideParsers{$offset}{$hdrName} = 1;
			}
			push (@$parserSeq, $hdrName, $offset);
			$pSeqStr .= " $hdrName " . $offset;
			$decLen = $currSeqLen + $hdrLastDecByte;
			$currSeqLen += $hdrLen;

			# Calculate the number of cycles to parse the current header
			my $decCyc = POSIX::ceil(($offset + $hdrLastDecByte) / $wordByteWidth);
			my $hdrCyc = POSIX::ceil(($offset + $hdrLen) / $wordByteWidth);
			$decCyc -= 1 if ($parseWidth > $wordByteWidth);
			$decCyc = 1 if ($hdrCyc > 0 && $decCyc == 0);

			if ($offset == 0) {
				$parseCycles += $decCyc;
			}

			if ($hdrCyc > $decCyc) {
				$parseCycles = max($parseCycles, POSIX::ceil(($hdrPos + $hdrLen) / $wordByteWidth));
			}
		}
		push (@parserSeqs, $parserSeq) if (!defined($pSeqStrs{$pSeqStr}));
		$pSeqStrs{$pSeqStr} = 1;

		# Stats tracking
		if ($lag > $maxLag) {
			$maxLag = $lag;
			$maxLagSeq = $hdrSeq;
		}
		my $hdrSeqLen = $hdrSeq->[-1]->{'pos'} + $hdrSeq->[-1]->{'len'};

		my $numPktRxBytes = max($hdrSeqLen, MIN_PKT_SIZE);
		my $numPktRxBytes_FCS_IPG = max($hdrSeqLen, MIN_PKT_SIZE);
		my $numPktRxCycles_no_FCS_IPG = POSIX::ceil($numPktRxBytes / $wordByteWidth);
		my $numPktRxCycles_with_FCS_IPG = POSIX::floor(($numPktRxBytes + FCS_LEN + IPG) / $wordByteWidth);
		my $numPktRxCycles = max($numPktRxCycles_no_FCS_IPG, $numPktRxCycles_with_FCS_IPG);
		my $pktParseRatio = $parseCycles / $numPktRxCycles;
		if ($pktParseRatio > $worstPktParseRatio) {
			$worstPktParseCycles = $parseCycles;
			$worstPktRxCycles = $numPktRxCycles;
			$worstPktParseSeq = $hdrSeq;
			$worstPktParseRatio = $pktParseRatio;
		}

		my $numHdrRxBytes = $hdrSeqLen;
		my $numHdrRxCycles = POSIX::ceil($numHdrRxBytes / $wordByteWidth);
		my $hdrParseRatio = $parseCycles / $numHdrRxCycles;
		if ($hdrParseRatio > $worstHdrParseRatio) {
			$worstHdrParseCycles = $parseCycles;
			$worstHdrRxCycles = $numHdrRxCycles;
			$worstHdrParseSeq = $hdrSeq;
			$worstHdrParseRatio = $hdrParseRatio;
		}

		if ($hdrSeqLen > $maxHdrSeqLen) {
			$maxHdrSeqLen = $hdrSeqLen;
			$maxHdrSeq = $hdrSeq;
		}
	}
	return \%parsers, \@parserSeqs, \%wideParsers, $maxLag, $maxLagSeq;
}

# Return the previously calculated parser locations
sub getParserLocs {
	return \%parsers, \@parserSeqs, \%wideParsers, $maxLag, $maxLagSeq;
}

# Return previously calculated worst-case packet parse info
sub getWorstPktInfo {
	return $worstPktParseCycles, $worstPktRxCycles, $worstPktParseSeq;
}

# Return previously calculated worst-case hdr parse info
sub getWorstHdrInfo {
	return $worstHdrParseCycles, $worstHdrRxCycles, $worstHdrParseSeq;
}

sub getMaxHdrSeqLen {
	return $maxHdrSeqLen;
}

sub getMaxHdrSeq {
	return $maxHdrSeq;
}

# Get the maximum number of parsers in a single cycle
sub maxParserSeqLen {
	my $maxParserSeqLen = 1;
	for (my $i = 0; $i < scalar(@parserSeqs); $i++) {
		my $parserSeq = $parserSeqs[$i];
		my $hdrName = $parserSeq->[0];
		my $offset = $parserSeq->[1];

		my $seqLen = scalar(@$parserSeq) / 2;
		$seqLen += 1 if $hdrName eq $firstHdrName and $firstHdrDecPos <= $wordByteWidth;
		$maxParserSeqLen = $seqLen if $seqLen > $maxParserSeqLen;
	}
	return $maxParserSeqLen;
}

# Record the set of headers in use
sub setHeaders {
	$headers = shift;
	$numHeaders = scalar keys %$headers;
}

# Get the list of headers
sub getHeaderNames {
	return sort(keys(%$headers));
}

# How many headers are there?
sub numHeaders {
	return $numHeaders;
}

# Get a particular header
sub getHeader {
	my $hdrName = shift;
	return $headers->{$hdrName};
}

# Get the first header that we expect
sub getFirstHdrName {
	return $firstHdrName;
}

# Get the first header that we expect
sub getFirstHdrDecPos {
	return $firstHdrDecPos;
}

# Get the first header that we expect
sub getFirstHdrEnum {
	my $hdr = getHeader($firstHdrName);
	return $hdr->{'enumName'};
}

# Compare two numbers
sub cmpNum ($$) {
	return $_[0] <=> $_[1];
}

# Sort a parser sequence
sub sortParserSeqLen {
	return sort parserSeqCompLen @_;
}

# Parser sequence comparison function
sub parserSeqCompLen {
	return scalar(@$b) <=> scalar(@$a) if scalar(@$b) != scalar(@$a);
	for (my $i = 0; $i < scalar(@$b); $i += 2) {
		my $aHdrName = $a->[$i];
		my $aOffset = $a->[$i + 1];
		my $bHdrName = $b->[$i];
		my $bOffset = $b->[$i + 1];
		my $ret = $aHdrName cmp $bHdrName;
		return $ret if $ret != 0;
		$ret = $aOffset <=> $bOffset;
		return $ret if $ret != 0;
	}
	return 0;
}

# Sort a parser sequence
sub sortParserSeqAlpha {
	return sort parserSeqCompAlpha @_;
}

# Parser sequence comparison function
sub parserSeqCompAlpha {
	for (my $i = 0; $i < scalar(@$a) and $i < scalar(@$b); $i += 2) {
		my $aHdrName = $a->[$i];
		my $aOffset = $a->[$i + 1];
		my $bHdrName = $b->[$i];
		my $bOffset = $b->[$i + 1];
		my $ret = $aHdrName cmp $bHdrName;
		return $ret if $ret != 0;
		$ret = $aOffset <=> $bOffset;
		return $ret if $ret != 0;
	}
	return scalar(@$a) <=> scalar(@$b);
}

# Build the parameters to pass into the generic processor
sub genericProcParams {
	my $hdr = shift;

	my %parameters;
        $parameters{'HdrLen'} = $hdr->{'len'};
        $parameters{'HdrLenBytes'} = $hdr->{'lenBytes'};
        my $hdrLenMap = [];
        foreach my $match (sort(keys(%{$hdr->{'lenMap'}}))) {
            push @$hdrLenMap, {'match' => $match, 'len' => $hdr->{'lenMap'}->{$match}};
        }
        $parameters{'HdrLenMap'} = $hdrLenMap;
        $parameters{'NxtHdrBytes'} = $hdr->{'nxtHdrBytes'};
        my $nxtHdrMap = [];
        foreach my $match (sort(keys(%{$hdr->{'nxtHdrMap'}}))) {
            my $nxtHdrEnum = "UNKNOWN";
            my $nxtHdrLen = 0;
            my $nxtHdr = getHeader($hdr->{'nxtHdrMap'}->{$match});
            if (defined($nxtHdr)) {
                $nxtHdrEnum = $nxtHdr->{'enumName'};
                $nxtHdrLen = $nxtHdr->{'len'};
            }
            push @$nxtHdrMap, {'match' => $match, 'nxt_hdr' => $nxtHdrEnum, 'nxt_hdr_len' => $nxtHdrLen};
        }
        $parameters{'NxtHdrMap'} = $nxtHdrMap;

	return %parameters;
}

# Generate a parser name
sub makeParserName {
	my ($hdrName, $offset) = @_;

	my $hdr = getHeader($hdrName);
	my $hdrEnum = $hdr->{'enumName'};

	return "${hdrEnum}_${offset}_parser";
}

# Generate a parser interface name
sub makeParserIfaceName {
	return makeParserName(@_) . "_ifc";
}

sub log2 {
    my $val = shift;

    return 1 if ($val <= 1) ;
    my $ret = 0;
    while ($val > 1) {
        $val /= 2;
        $ret++;
    }
    return $ret;
}

sub log10 {
    my $val = shift;

    return 1 if ($val <= 10) ;
    my $ret = 0;
    while ($val > 1) {
        $val /= 10;
        $ret++;
    }
    return $ret;
}

sub findHashMatch {
	my ($hash, $val) = @_;

	my $retKey;
	foreach my $key (sort(keys(%$hash))) {
		my $hashVal = $hash->{$key};
		$retKey = $key if !defined($retKey) and $val eq $hashVal;
	}
	return $retKey;
}

sub extractMatchBytes {
	my $matchVal = shift;

	# Attempt to split the match value
	$matchVal =~ /^(\d+)'b(.*)$/;
	my @matchBytesStr = ( $2 =~ m/......../g);
	my @matchBytes;
	foreach my $byteStr (@matchBytesStr) {
		my $byte = 0;
		for (my $i = 0; $i < 8; $i++) {
			$byte += 2 ** (7 - $i) * (substr($byteStr, $i, 1) eq '1' ? 1 : 0);
		}
		push @matchBytes, $byte;
	}

	return @matchBytes;
}

sub extractMatchMaskBytes {
	my $matchVal = shift;

	# Attempt to split the match value
	$matchVal =~ /^(\d+)'b(.*)$/;
	my @matchBytesStr = ( $2 =~ m/......../g);
	my @matchBytes;
	my @maskBytes;
	foreach my $byteStr (@matchBytesStr) {
		my $match = 0;
		my $mask = 0;
		for (my $i = 0; $i < 8; $i++) {
			$match += 2 ** (7 - $i) * (substr($byteStr, $i, 1) eq '1' ? 1 : 0);
			$mask += 2 ** (7 - $i) * (substr($byteStr, $i, 1) ne '?' ? 1 : 0);
		}
		push @matchBytes, $match;
		push @maskBytes, $mask;
	}

	return (\@matchBytes, \@maskBytes);
}

sub genPktByteStream {
	my $hdrSeq = shift;
	my $randomData = shift;
	my $seqData = shift;
	my @bytes;
	my @mask;

	# Walk through the headers, adding bytes to the byte array
	my @remainder = ();
	my @remainderMask = ();
	for (my $i = 0; $i < scalar(@$hdrSeq); $i++) {
		my $hdrName = $hdrSeq->[$i]->{'name'};
		my $hdrPos = $hdrSeq->[$i]->{'pos'};
		my $hdrLen = $hdrSeq->[$i]->{'len'};
		my $hdrLastDecByte = getLastDecByte($hdrName) + 1;
		my $byteCnt = max($hdrLen, $hdrLastDecByte, scalar(@remainder));
		my $hdr = $headers->{$hdrName};

		my $nxtHdrName;
		#my $nxtHdrPos;
		#my $nxtHdrLen;
		if ($i < scalar(@$hdrSeq) - 1) {
			$nxtHdrName = $hdrSeq->[$i + 1]->{'name'};
		}

		my @hdrBytes = (0) x $byteCnt;
		my @hdrMask = (255) x $byteCnt;

		# Copy across any remainder
		for (my $j = 0; $j < scalar(@remainder); $j++) {
			$hdrBytes[$j] = $remainder[$j];
			$hdrMask[$j] = $remainderMask[$j];
		}

		# Populate the length
		if ($hdr->{'len'} == 0 && scalar(keys(%{$hdr->{'lenMap'}})) > 0) {
			my $matchVal = findHashMatch($hdr->{'lenMap'}, $hdrLen);
			confess "Cannot find length '$hdrLen' in header '$hdrName'" unless defined($matchVal);
			my @matchBytes = extractMatchBytes($matchVal);

			my $lenBytes = $hdr->{'lenBytes'};
			for (my $j = 0; $j < scalar(@$lenBytes); $j++) {
				$hdrBytes[$lenBytes->[$j]] |= $matchBytes[$j];
				$hdrMask[$lenBytes->[$j]] = 0;
			}
		}

		# Populate the match
		my @matchBytes;
		if (defined($nxtHdrName)) {
			my $matchVal = findHashMatch($hdr->{'nxtHdrMap'}, $nxtHdrName);
			confess "Cannot find next header '$nxtHdrName' in header '$hdrName'" unless defined($matchVal);
			@matchBytes = extractMatchBytes($matchVal);

			my $nxtHdrBytes = $hdr->{'nxtHdrBytes'};
			for (my $j = 0; $j < scalar(@$nxtHdrBytes); $j++) {
				$hdrBytes[$nxtHdrBytes->[$j]] |= $matchBytes[$j];
				$hdrMask[$nxtHdrBytes->[$j]] = 0;
			}
		} elsif (defined($headers->{$hdrName}->{'defNxtHdrVal'})) {
			@matchBytes = extractMatchBytes($headers->{$hdrName}->{'defNxtHdrVal'});

		} else {
			@matchBytes = (0) x scalar(@{$hdr->{'nxtHdrBytes'}});
		}
		my $nxtHdrBytes = $hdr->{'nxtHdrBytes'};
		for (my $j = 0; $j < scalar(@$nxtHdrBytes); $j++) {
			$hdrBytes[$nxtHdrBytes->[$j]] |= $matchBytes[$j];
			$hdrMask[$nxtHdrBytes->[$j]] = 0;
		}

		# Add the header to the byte sequence
		push @bytes, @hdrBytes[0 .. $hdrLen-1];
		push @mask, @hdrMask[0 .. $hdrLen-1];

		# Store any leftovers in the remainder array
		if (scalar(@hdrBytes) > $hdrLen) {
			@remainder = @hdrBytes[$hdrLen .. scalar(@hdrBytes) - 1];
			@remainderMask = @hdrMask[$hdrLen .. scalar(@hdrBytes) - 1];
		} else {
			@remainder = ();
			@remainderMask = ();
		}
	}

	# Append any remainder
	push @bytes, @remainder;
	push @mask, @remainderMask;

	# Pad the byte sequence if needed
	if (scalar(@bytes) < MIN_PKT_SIZE) {
		my $totalHdrLen = scalar(@bytes);
		my @payload = (0) x (MIN_PKT_SIZE - $totalHdrLen);
		my @payloadMask = (255) x (MIN_PKT_SIZE - $totalHdrLen);
		push @bytes, @payload;
		push @mask, @payloadMask;
	}

	if ($randomData || $seqData) {
		my @bgData = (0 .. scalar(@bytes));
		if ($randomData) {
			@bgData = map(POSIX::floor(rand(256)), @bgData);
		} else {
			@bgData = map { $_ % 256; } @bgData;
		}
		for (my $i = 0; $i < scalar(@bytes); $i++) {
			$bytes[$i] = $bytes[$i] + ($bgData[$i] & $mask[$i]);
		}
	}

	return @bytes;
}

# Generate the extracted field byte stream
sub genExtractByteStream {
	my $hdrSeq = shift;
	my $pkt = shift;
	my $evByteWidth = shift;

	#my %offsets;
	my @extFields = (0) x $evByteWidth;
	my @extFieldsVld = (0) x $evByteWidth;

	foreach my $hdrInfo (@$hdrSeq) {
		my $hdrName = $hdrInfo->{'name'};
		my $hdrPos = $hdrInfo->{'pos'};
		my $hdrLen = $hdrInfo->{'len'};
		my $hdrMap = $extractMap{$hdrName};
		my $pktPos = 0;
		my $offsetName = $hdrName;


		#$offsets{$hdrName} = 0 if !defined($offsets{$hdrName});
		#my $offset = $offsets{$hdrName};
		my $offset = ($hdrInfo->{'cnt'} - 1) * $extractByteCnt{$hdrName};
		#$offsets{$hdrName} += $extractByteCnt{$hdrName};

		foreach my $extractRow (@$hdrMap) {
			foreach my $extractDst (@$extractRow) {
				if ($extractDst != 0) {
					$extFields[$extractDst-1 + $offset] = $pkt->[$pktPos + $hdrPos];
					$extFieldsVld[$extractDst-1 + $offset] = 1;
				}
				$pktPos++;
			}
		}
	}

	return \@extFields, \@extFieldsVld;
}

# Print a verilog array
sub printVerilogArray {
	my $arr = shift;
	my $name = shift;
	my $fmtStr = shift;
	my $indent = shift;
	my $unpackedInit = shift;

	$indent = 0 if !defined($indent);
	$unpackedInit = 0 if !defined($unpackedInit);
	my $arrLen = scalar(@$arr);

	print " " x $indent .  "assign $name = ";
	print "'" if ($unpackedInit);
	print "{\n";
	for (my $j = 0; $j < $arrLen; $j++) {
		print " " x $indent . "   " if ($j % 8 == 0);
		printf " $fmtStr", $arr->[$j];
		print "," if ($j != $arrLen - 1);
		if ($j % 8 == 7) {
			print "\n";
		}
	}
	print "\n" if ($arrLen % 8 != 0);
	print ' ' x $indent . "};\n";
}

# Store the list of header sequences
sub setHdrSeqs {
	$hdrSeqs = shift;
}

# Return the header sequences
sub getHdrSeqs {
	return $hdrSeqs;
}

sub minPktSize {
	return MIN_PKT_SIZE;
}

sub genHdrPosCheck {
	my ($dblWide, $pos) = @_;
	my $result = "";
	if ($dblWide && $pos < 2 * $wordByteWidth) {
		$result .= "((ifc.hdr_pos == " . wordStartPos($pos) . ") ||";
		$result .= " (ifc.hdr_pos == '0 && ifc.this_dbl_wide))";
	} else {
		$result .= "ifc.hdr_pos == " . wordStartPos($pos);
	}

	return $result;
}

sub calcOffsetsFromLen {
    my $matchLen = shift;
    my $lastExtractByte = shift;
    my $offset = shift;
    my $dblWide = shift;

    my $matchLenAdj = $matchLen + $offset;
    my $lastExtractByteAdj = $lastExtractByte + $offset;

    my $skipBytes = wordStartPos($matchLenAdj - 1) != wordStartPos($lastExtractByteAdj) ? 1 : 0;
    my $skipBytesDW = ($lastExtractByteAdj + 1) <= 2 * $wordByteWidth &&
        $matchLenAdj > 2 * $wordByteWidth;
    my $offsetInc = -$matchLenAdj % $wordByteWidth;
    if ($skipBytes) {
        $offsetInc = $lastExtractByteAdj + 1 + $wordByteWidth - 1;
        $offsetInc -= $offsetInc % $wordByteWidth + $matchLenAdj;
    }
    my $offsetIncDW = -$matchLenAdj % (2 * $wordByteWidth);
    if ($skipBytesDW) {
        $offsetIncDW = $lastExtractByteAdj + 1 + (2 * $wordByteWidth) - 1;
        $offsetIncDW -= $offsetIncDW % (2 * $wordByteWidth) + $matchLenAdj;
    }

    # Need to jump from start of last extract word to the end of the header
    my $posInc = $matchLenAdj;
    my $posIncDW = $matchLenAdj;

    # A header decision point will either fit within the word, or will start at offset=0
    if ($offset == 0) {
	    my $firstPosInLastDecWord = int($lastExtractByte / $wordByteWidth) * $wordByteWidth;
	    $posInc -= $firstPosInLastDecWord;
	    # posIncDW is only used when the decision word is the first word,
	    # so no need to subtract the firstPosInLastDecWord
    }

    # Clear posIncDW if the decision point stretches beyond the end of the current double-word
    # (Only the first word is processed in double-wide mode.)
    if ($lastExtractByte + $offset > 2 * $wordByteWidth) {
	$posIncDW = 0;
    }

    # decRem = # of unprocessed bytes in decision word
    # +ve value means that unprocessed bytes remain
    # -ve value means that we are skipping bytes
    my $extractLen = POSIX::floor(($lastExtractByteAdj + $wordByteWidth) / $wordByteWidth) * $wordByteWidth;
    my $decRem = $extractLen - $matchLenAdj;
    my $extractLenDW = $extractLen;
    if ($extractLenDW < 2 * $wordByteWidth) {
	    $extractLenDW = 2 * $wordByteWidth;
    }
    my $decRemDW = $extractLenDW - $matchLenAdj;

    my $shiftAmtInc = $decRem % $wordByteWidth;
    my $bufRdAmt = -($decRem - $shiftAmtInc) / $wordByteWidth;

    my ($shiftAmtIncDW, $bufRdAmtDW);
    if ($decRemDW >= 0) {
	    $shiftAmtIncDW = $decRemDW % (2 * $wordByteWidth);
	    $bufRdAmtDW = -($decRemDW - $shiftAmtInc) / (2 * $wordByteWidth);
    } else {
	    $shiftAmtIncDW = $decRemDW % $wordByteWidth;
	    $bufRdAmtDW = -($decRemDW - $shiftAmtInc) / $wordByteWidth;
    }

    # Increment read amount by one (read the current word by default)
    $bufRdAmt += 1;
    $bufRdAmtDW += 1;

    if ($dblWide) {
	    if ($bufRdAmt == 0) {
		    $bufRdAmt = 1;
		    $shiftAmtInc += $wordByteWidth;
	    }
	    if ($bufRdAmtDW == 0) {
		    $bufRdAmtDW = 1;
		    $shiftAmtIncDW += $wordByteWidth;
	    }
    }

    return ($offsetInc, $offsetIncDW, $posInc, $posIncDW, $shiftAmtInc, $shiftAmtIncDW, $bufRdAmt, $bufRdAmtDW);
}

sub printLenCaseBranch {
    my $matchStr = shift;
    my $matchLen = shift;
    my $lastExtractByte = shift;
    my $offset = shift;
    my $dblWideEn = shift;
    my $thisDblWide = shift;

    my ($offsetInc, $offsetIncDW, $posInc, $posIncDW, $shiftAmtInc, $shiftAmtIncDW, $bufRdAmt, $bufRdAmtDW) = calcOffsetsFromLen($matchLen, $lastExtractByte, $offset, $dblWideEn);

    my $lastExtractByteAdj = $lastExtractByte + $offset;
    my $hdrPosCheck = genHdrPosCheck($thisDblWide, $lastExtractByteAdj);
    my $lastWordLen = $matchLen % $wordByteWidth;
    $lastWordLen = $wordByteWidth if $wordByteWidth == 0;
    my $lastWordPos = POSIX::floor(($matchLen - 1) / $wordByteWidth) * $wordByteWidth;

    $matchLen = fmtWidthVal($matchLen);
    $lastWordLen = fmtWidthVal($lastWordLen);
    $lastWordPos = fmtWidthVal($lastWordPos);
    $posInc = fmtWidthVal($posInc);
    $posIncDW = fmtWidthVal($posIncDW);
    $offsetInc = fmtWidthVal($offsetInc);
    $offsetIncDW = fmtWidthVal($offsetIncDW);
    $bufRdAmt = fmtBRAmt($bufRdAmt);
    $bufRdAmtDW = fmtBRAmt($bufRdAmtDW);
    $shiftAmtInc = fmtVal($shiftAmtInc, getShiftWidth());
    $shiftAmtIncDW = fmtVal($shiftAmtIncDW, getShiftWidth());

    print <<__LEN_CASE_BRANCH__;

        $matchStr : begin
            ifc.curr_hdr_len = $matchLen;
            ifc.curr_hdr_last_word_len = $lastWordLen;
            ifc.curr_hdr_last_word_pos = $lastWordPos;
            ifc.last_word = ifc.data_vld && $hdrPosCheck;
            ifc.nxt_hdr_proc_pos_inc = ifc.all_dbl_wide ? $posIncDW : $posInc;
            ifc.nxt_hdr_offset_inc = ifc.all_dbl_wide ? $offsetIncDW : $offsetInc;
            ifc.buf_rd_amt = ifc.all_dbl_wide ? $bufRdAmtDW : $bufRdAmt;
            ifc.shift_amt_inc = ifc.all_dbl_wide ? $shiftAmtIncDW : $shiftAmtInc;
        end
__LEN_CASE_BRANCH__

}

# Get the minimum header length
sub getMinHdrLen {
	my $hdr = shift;
	$hdr = $headers->{$hdr} unless ref($hdr);

	my $minHdrLen = $hdr->{'len'};
	if ($minHdrLen == 0 && scalar(keys(%{$hdr->{'lenMap'}})) > 0) {
		$minHdrLen = min(values(%{$hdr->{'lenMap'}}));
	}

	return $minHdrLen;
}

# Calculate the tables for extracting fields
sub calcExtractTables {
	# Assume that the headers have already been set via setHeaders

	# Reset vars
	%extractByteCnt = ();
	%extractMap = ();
	%extractHdrCnt = ();
	$totalExtractBytes = 0;
	$totalExtractEntries = 1; # Start at 1 because of DONE state

	# Walk through the list of all header sequences and count the numbers of each header
	my %hdrSeen;
	my @hdrOrder;
	foreach my $hdrSeq (@$hdrSeqs) {
		my %currCnt;

		# Walk through the headers in a single header list and count header occurences
		foreach my $hdrInfo (@$hdrSeq) {
			my $hdrName = $hdrInfo->{'name'};

			if (!defined($currCnt{$hdrName})) {
				$currCnt{$hdrName} = 0;
			}
			$currCnt{$hdrName}++;

			# Update the hdr order list
			if (!defined($hdrSeen{$hdrName})) {
				$hdrSeen{$hdrName} = 1;
				push @hdrOrder, $hdrName;
			}
		}

		# Update the header counts
		foreach my $hdrName (keys(%currCnt)) {
			if (!defined($extractHdrCnt{$hdrName}) or $currCnt{$hdrName} > $extractHdrCnt{$hdrName}) {
				$extractHdrCnt{$hdrName} = $currCnt{$hdrName};
			}
		}
	}

	# Create a sorted list of headers. Sort order is as follows:
	#  1. first header
	#  2. children of first header
	#  3. children of second headers
	#  4. etc
	my $firstHdr;
	my %succ;
	foreach my $hdrSeq (@$hdrSeqs) {
		my $prevHdr;
		foreach my $hdrInfo (@$hdrSeq) {
			my $hdrName = $hdrInfo->{'name'};
			
			$firstHdr = $hdrName if !defined($firstHdr);

			if (defined($prevHdr)) {
				$succ{$prevHdr} = {} if !defined($succ{$prevHdr});
				$succ{$prevHdr}->{$hdrName} = 1;
			}
			$prevHdr = $hdrName;
		}
	}

	@hdrOrder = ();
	my %seen;
	my @pending;
	push @pending, $firstHdr;
	$seen{$firstHdr} = 1;
	while (scalar(@pending) > 0) {
		my $hdrName = shift(@pending);
		push @hdrOrder, $hdrName;

		if (defined($succ{$hdrName})) {
			my @nxtHdrs = sort(keys(%{$succ{$hdrName}}));

			foreach my $hdr (@nxtHdrs) {
				if (!defined($seen{$hdr})) {
					push @pending, $hdr;
					$seen{$hdr} = 1;
				}
			}
		}
	}

	# Walk through the headers and work out how many bytes are required for each
	my $target = 1;
	foreach my $hdrName (@hdrOrder) {
		my $hdr = $headers->{$hdrName};
		my $minHdrLen = getMinHdrLen($hdr);

		# Identify the bytes needing extraction
		my @extractBytes = (0) x $minHdrLen;
		my $lastExtractByte = -1;
		foreach my $fieldInfo (@{$hdr->{'extractFields'}}) {
			my $fn = $fieldInfo->{'name'};
			my $fp = $fieldInfo->{'pos'};
			my $fl = $fieldInfo->{'len'};

			$fn = "$hdrName:$fn";
			for (my $i = int($fp / 8); $i <= int(($fp + $fl - 1) / 8); $i++) {
				$extractBytes[$i] = 1;
				$lastExtractByte = $i if ($i > $lastExtractByte);
			}
		}

		# Create an extraction map
		# (ie. map of bytes to extract and locations to send them to)
		my $cnt = 0;
		my $extractByteMap = [];
		for (my $i = 0; $i <= $lastExtractByte || $i == 0; $i += $wordByteWidth) {
			my $nullMap = [(0) x $wordByteWidth];
			push @$extractByteMap, $nullMap;
			$totalExtractEntries++;
		}
		for (my $i = 0; $i < $minHdrLen; $i++) {
			if ($extractBytes[$i]) {
				$extractByteMap->[int($i / $wordByteWidth)]->[$i % $wordByteWidth] = $target++;
				$cnt++;
			}
		}
		$extractByteCnt{$hdrName} = $cnt;
		$target += $cnt * ($extractHdrCnt{$hdrName} - 1);
		$extractMap{$hdrName} = $extractByteMap;
	}
	$totalExtractBytes = $target;

	# Add dummy entries for unused headers
	foreach my $hdrName (keys(%$headers)) {
		if (!defined($extractMap{$hdrName})) {
			$extractByteCnt{$hdrName} = 0;
			$extractMap{$hdrName}[0] = [(0) x $wordByteWidth];
			$extractHdrCnt{$hdrName} = 0;
			$totalExtractEntries++;
		}
	}

	return ($totalExtractBytes, $totalExtractEntries, \%extractMap, \%extractByteCnt, \%extractHdrCnt);
}

# Return the already calculated extract tables
sub getExtractTables {
	return ($totalExtractBytes, $totalExtractEntries, \%extractMap, \%extractByteCnt, \%extractHdrCnt);
}

sub calcIfcArrays {
	# Start by creating the maps
	my @offsets = sort(keys(%parsers));
	foreach my $offset (@offsets) {
		foreach my $hdrName (sort(keys(%{$parsers{$offset}}))) {
			my $hdrDblWide = defined($wideParsers{$offset}->{$hdrName});
			$parserIfcMap{$offset} = {} if !defined($parserIfcMap{$offset});
			if ($hdrDblWide) {
				$parserIfcMap{$offset}->{$hdrName} = $numIfcDW++;
			} else {
				$parserIfcMap{$offset}->{$hdrName} = $numIfcSW++;
			}
		}
	}
}

sub getIfcArrayPos {
	my ($offset, $hdrName) = @_;
	return $parserIfcMap{$offset}->{$hdrName};
}

sub createIfcInsts {
	my ($ifc, $ifcDblWide) = @_;

	my @offsets = sort(keys(%parsers));
	foreach my $offset (@offsets) {
		foreach my $hdrName (sort(keys(%{$parsers{$offset}}))) {
			my $hdr = getHeader($hdrName);
			my $hdrEnum = $hdr->{'enumName'};
			my $ifc_name = "${hdrEnum}_${offset}_parser_ifc";
			my $newIfc;
			if (defined($wideParsers{$offset}->{$hdrName})) {
				$newIfc = $ifcDblWide->clone($ifc_name);
			} else {
				$newIfc = $ifc->clone($ifc_name);
			}
			print $newIfc->instantiate() . "( .* );\n";
		}
	}
}

sub ifcArrayCopy {
	my $swap = shift;

	my @offsets = sort(keys(%parsers));
	foreach my $offset (@offsets) {
		foreach my $hdrName (sort(keys(%{$parsers{$offset}}))) {
			my $hdr = getHeader($hdrName);
			my $hdrEnum = $hdr->{'enumName'};
			my $dstIfc = "${hdrEnum}_${offset}_parser_ifc";
			my $srcIfc;
			if (defined($wideParsers{$offset}->{$hdrName})) {
				$srcIfc = "parser_ifcs_wide[" . getIfcArrayPos($offset, $hdrName) . "]";
			} else {
				$srcIfc = "parser_ifcs[" . getIfcArrayPos($offset, $hdrName) . "]";
			}
			if ($swap) {
				my $tmp = $srcIfc;
				$srcIfc = $dstIfc;
				$dstIfc = $tmp;
			}
			print <<__ARRAY_COPY__;
assign $dstIfc.nxt_hdr_info_vld       = $srcIfc.nxt_hdr_info_vld;
assign $dstIfc.nxt_hdr_type           = $srcIfc.nxt_hdr_type;
assign $dstIfc.nxt_hdr_len            = $srcIfc.nxt_hdr_len;
assign $dstIfc.nxt_hdr_proc_pos_inc   = $srcIfc.nxt_hdr_proc_pos_inc;
assign $dstIfc.nxt_hdr_offset_inc     = $srcIfc.nxt_hdr_offset_inc;
assign $dstIfc.curr_hdr_info_vld      = $srcIfc.curr_hdr_info_vld;
assign $dstIfc.curr_hdr_len           = $srcIfc.curr_hdr_len;
assign $dstIfc.curr_hdr_last_word_len = $srcIfc.curr_hdr_last_word_len;
assign $dstIfc.curr_hdr_last_word_pos = $srcIfc.curr_hdr_last_word_pos;
assign $dstIfc.last_word              = $srcIfc.last_word;
assign $dstIfc.buf_rd_amt             = $srcIfc.buf_rd_amt;
assign $dstIfc.shift_amt_inc          = $srcIfc.shift_amt_inc;

__ARRAY_COPY__
		}
	}
}

sub copyIfcToArray {
	ifcArrayCopy(1);
}

sub copyIfcFromArray {
	ifcArrayCopy(0);
}

sub getIfcArrayLen {
	return $numIfcSW;
}

sub getIfcArrayLenDW {
	return $numIfcDW;
}

sub getFinalParsers {
	my $parsersSW = {};
	my $nonFinalParsersSW = {};
	my $parsersDW = {};
	my $nonFinalParsersDW = {};

	foreach my $parserSeq (@parserSeqs) {
		my $prevHdrNameSW;
		my $prevOffsetSW;

		my $prevHdrNameDW;
		my $prevOffsetDW;

		for (my $i = 0; $i < scalar(@$parserSeq); $i += 2) {
			my $hdrName = $parserSeq->[$i];
			my $offset = $parserSeq->[$i + 1];
			my $lastHdrDecPos = getLastDecByte($hdrName);

			if ($i > 0) {
				$nonFinalParsersSW->{$prevOffsetSW} = {} if !defined($nonFinalParsersSW->{$prevOffsetSW});
				$nonFinalParsersSW->{$prevOffsetSW}->{$prevHdrNameSW} = 1;

				$nonFinalParsersDW->{$prevOffsetDW} = {} if !defined($nonFinalParsersDW->{$prevOffsetDW});
				$nonFinalParsersDW->{$prevOffsetDW}->{$prevHdrNameDW} = 1;
			}

			if ($offset == 0 || $offset + $lastHdrDecPos < $wordByteWidth) {
				$parsersSW->{$offset} = {} if !defined($parsersSW->{$offset});
				$parsersSW->{$offset}->{$hdrName} = 1;

				$prevHdrNameSW = $hdrName;
				$prevOffsetSW = $offset;
			}

			$parsersDW->{$offset} = {} if !defined($parsersDW->{$offset});
			$parsersDW->{$offset}->{$hdrName} = 1;

			$prevHdrNameDW = $hdrName;
			$prevOffsetDW = $offset;
		}
	}

	# Delete any non-final parsers
	foreach my $offset (keys(%$nonFinalParsersDW)) {
		foreach my $hdrName (keys(%{$nonFinalParsersDW->{$offset}})) {
			if (defined($nonFinalParsersSW->{$offset}) and
			    defined($nonFinalParsersSW->{$offset}->{$hdrName}) and
			    defined($parsersSW->{$offset}) and
			    defined($parsersSW->{$offset}->{$hdrName})) {
				delete $parsersSW->{$offset}->{$hdrName};
			}
			if (defined($parsersDW->{$offset}) and
			    defined($parsersDW->{$offset}->{$hdrName})) {
				delete $parsersDW->{$offset}->{$hdrName};
			}
		}
	}

	# Delete any empty offsets
	foreach my $offset (keys(%$parsersDW)) {
		if (scalar(keys(%{$parsersSW->{$offset}})) == 0) {
			delete $parsersSW->{$offset};
		}
		if (scalar(keys(%{$parsersDW->{$offset}})) == 0) {
			delete $parsersDW->{$offset};
		}
	}

	return $parsersSW, $parsersDW;
}

sub getPrecedingParsers {
	my $precedingParsers = {};

	foreach my $parserSeq (@parserSeqs) {
		my $prevHdrName = $parserSeq->[0];
		my $prevOffset = $parserSeq->[1];

		for (my $i = 2; $i < scalar(@$parserSeq); $i += 2) {
			my $hdrName = $parserSeq->[$i];
			my $offset = $parserSeq->[$i + 1];

			my $curr = $precedingParsers;

			$curr->{$offset} = {} if !defined($curr->{$offset});
			$curr = $curr->{$offset};

			$curr->{$hdrName} = {} if !defined($curr->{$hdrName});;
			$curr = $curr->{$hdrName};

			$curr->{$prevOffset} = {} if !defined($curr->{$prevOffset});;
			$curr = $curr->{$prevOffset};

			$curr = $curr->{$prevHdrName} = 1;

			$prevHdrName = $hdrName;
			$prevOffset = $offset;
		}
	}

	return $precedingParsers;
}

sub getMaxHdrLen {
	return MAX_HDR_LEN;
}

sub getHdrWidth {
	return log2(getMaxHdrLen());
}

sub fmtWidthVal {
	my $val = shift;

	return fmtVal($val, getHdrWidth());
}

sub fmtVal {
	my $val = shift;
	my $width = shift;
	my $sign = $val >= 0 ? '' : '-';
	$val = abs($val);

	return sprintf("%s%d'd%d", $sign, $width, $val);
}

sub getBufRdWidth {
	return log2(POSIX::ceil(getMaxHdrLen() / $wordByteWidth)) + 1;
}

sub fmtBRAmt {
	my $val = shift;

	return fmtVal($val, getBufRdWidth());
}

sub setShiftWidth {
	$shiftWidth = shift;
}

sub getShiftWidth {
	return $shiftWidth;
}

# Calculate the tables for extracting fields
sub calcProgLookupTable {
	# Assume that the headers have already been set via setHeaders

	# Reset vars
	%progTableEntries = ();
	$progNumEntries = 0;
	$progUniqueEntries = 0;
	$progFirstLookupWords = undef;

	# Local vars
	my %progExtractByteCnt = ();
	my %progExtractMap = ();
	my %progExtractHdrCnt = ();
	my %progExtractHdrBase = ();
	my $progTotalExtractBytes = 0;
	my %progLookupBytes = ();

	my %progMaskSeqs = ();
	my %progMatchSeqs = ();
	my %progLens = ();
	my %progNxtHdrs = ();

	# Walk through the list of all header sequences and count the numbers of each header
	foreach my $hdrSeq (@$hdrSeqs) {
		my %currCnt;

		# Walk through the headers in a single header list and count header occurences
		foreach my $hdrInfo (@$hdrSeq) {
			my $hdrName = $hdrInfo->{'name'};

			if (!defined($currCnt{$hdrName})) {
				$currCnt{$hdrName} = 0;
			}
			$currCnt{$hdrName}++;
		}

		# Update the header counts
		foreach my $hdrName (keys(%currCnt)) {
			if (!defined($progExtractHdrCnt{$hdrName}) or $currCnt{$hdrName} > $progExtractHdrCnt{$hdrName}) {
				$progExtractHdrCnt{$hdrName} = $currCnt{$hdrName};
			}
		}
	}

	# Create a sorted list of headers. Sort order is as follows:
	#  1. first header
	#  2. children of first header
	#  3. children of second headers
	#  4. etc
	my $firstHdr;
	my %succ;
	foreach my $hdrSeq (@$hdrSeqs) {
		my $prevHdr;
		foreach my $hdrInfo (@$hdrSeq) {
			my $hdrName = $hdrInfo->{'name'};
			
			$firstHdr = $hdrName if !defined($firstHdr);

			if (defined($prevHdr)) {
				$succ{$prevHdr} = {} if !defined($succ{$prevHdr});
				$succ{$prevHdr}->{$hdrName} = 1;
			}
			$prevHdr = $hdrName;
		}
	}

	my @hdrOrder;
	my %seen;
	my @pending;
	push @pending, $firstHdr;
	$seen{$firstHdr} = 1;
	while (scalar(@pending) > 0) {
		my $hdrName = shift(@pending);
		push @hdrOrder, $hdrName;

		if (defined($succ{$hdrName})) {
			my @nxtHdrs = sort(keys(%{$succ{$hdrName}}));

			foreach my $hdr (@nxtHdrs) {
				if (!defined($seen{$hdr})) {
					push @pending, $hdr;
					$seen{$hdr} = 1;
				}
			}
		}
	}

	# Walk through the headers and work out how many bytes are required for
	# each and which bytes are used in the lookup
	my $target = 0;
	foreach my $hdrName (@hdrOrder) {
		my $hdr = $headers->{$hdrName};
		my $minHdrLen = getMinHdrLen($hdr);

		# Identify the bytes needing extraction
		my @extractBytes = (0) x $minHdrLen;
		my $lastExtractByte = -1;
		my $numExtractBytes = 0;
		foreach my $fieldInfo (@{$hdr->{'extractFields'}}) {
			my $fn = $fieldInfo->{'name'};
			my $fp = $fieldInfo->{'pos'};
			my $fl = $fieldInfo->{'len'};

			$fn = "$hdrName:$fn";
			for (my $i = int($fp / 8); $i <= int(($fp + $fl - 1) / 8); $i++) {
				$numExtractBytes += 1 if $extractBytes[$i] == 0;
				$extractBytes[$i] = 1;
				$lastExtractByte = $i if ($i > $lastExtractByte);
			}
		}

		# Create an extraction map
		# (ie. map of bytes to extract and locations to send them to)
		my $cnt = 0;
		my $extractByteMap = {};
		for (my $i = 0; $i < $minHdrLen; $i++) {
			if ($extractBytes[$i]) {
				$extractByteMap->{$i} = $cnt++;
			}
		}

		# Create a sorted list of the bytes used in lookups
		my %lookupByteMap;
		my $lenBytes = $hdr->{'lenBytes'};
		for (my $i = 0; $i < scalar(@$lenBytes); $i++) {
			$lookupByteMap{int($lenBytes->[$i])} = 1;
		}
		my $nxtHdrBytes = $hdr->{'nxtHdrBytes'};
		for (my $i = 0; $i < scalar(@$nxtHdrBytes); $i++) {
			$lookupByteMap{int($nxtHdrBytes->[$i])} = 1;
		}
		my @lookupBytes = sort(keys(%lookupByteMap));

		$progExtractByteCnt{$hdrName} = $cnt;
		$progExtractHdrBase{$hdrName} = $target;
		$target += $cnt * $extractHdrCnt{$hdrName};
		$progExtractMap{$hdrName} = $extractByteMap;
		$progLookupBytes{$hdrName} = \@lookupBytes;
	}
	$progTotalExtractBytes = $target;

	# Add dummy entries for unused headers
	foreach my $hdrName (keys(%$headers)) {
		if (!defined($extractMap{$hdrName})) {
			$progExtractByteCnt{$hdrName} = 0;
			$progExtractMap{$hdrName} = {};
			$progExtractHdrCnt{$hdrName} = 0;
			$progExtractHdrBase{$hdrName} = $target;
		}
	}

	my %progLGroups = ();
	my %progESGroups = ();
	my %progEDGroups = ();
	my %progShifts = ();

	# Work out the groups/shifts for each header type
	# (Naively partition into maximum-sized groups)
	foreach my $hdrName (@hdrOrder) {
		my $hdr = $headers->{$hdrName};
		my $minHdrLen = getMinHdrLen($hdr);

		my $lgroups = [];
		my $esgroups = [];
		my $edgroups = [];
		my $shifts = [];

		my @extractBytes = sort(cmpNum keys(%{$progExtractMap{$hdrName}}));
		my @lookupBytes = @{$progLookupBytes{$hdrName}};
		my $numExtractBytes = scalar(@extractBytes);
		my $numLookupBytes = scalar(@lookupBytes);
		my $ep = 0;
		my $lp = 0;
		my $startPos = 0;
		my $lgroup = [];
		my $esgroup = [];
		my $edgroup = [];
		my $lsize = 0;
		my $esize = 0;
		while ($ep < $numExtractBytes or $lp < $numLookupBytes) {
			my $nextEB = MAX_HDR_LEN + 1;
			my $nextLB = MAX_HDR_LEN + 1;
			$nextEB = $extractBytes[$ep] if ($ep < $numExtractBytes);
			$nextLB = $lookupBytes[$lp] if ($lp < $numLookupBytes);

			my $nextB = min($nextEB, $nextLB);

			# Record the groups
			if ($nextB >= $startPos + $progMaxRange or
				($nextB == $nextEB and $esize >= $progExtractWidth) or
				($nextB == $nextLB and $lsize >= $progLookupInputs)) {
				#print STDERR "PUSH\n";
				my $newStartPos = min($nextB, $startPos + $progMaxRange);

				push @$lgroups, $lgroup;
				push @$esgroups, $esgroup;
				push @$edgroups, $edgroup;
				push @$shifts, $newStartPos - $startPos;

				$startPos = $newStartPos;
				$lgroup = [];
				$esgroup = [];
				$edgroup = [];
				$lsize = 0;
				$esize = 0;

				# Jump to the next loop iteration in case of large breaks in headers
				next;
			}

			# Extract byte comes first
			if ($nextLB <= $nextEB) {
				if ($nextLB - $startPos + $progLookupWidth > $progMaxRange) {
					$nextLB = $startPos + $progLookupWidth - $progMaxRange;
				}
				push @$lgroup, $nextLB - $startPos;
				$lsize++;

				while ($lp < $numLookupBytes and $lookupBytes[$lp] < $nextLB + $progLookupWidth) {
					$lp++;
				}
			}
			else {
				# $nextEB < $nextLB
				push @$esgroup, $nextEB - $startPos;
				push @$edgroup, $progExtractHdrBase{$hdrName} + $progExtractMap{$hdrName}->{$nextEB};
				$esize++;
				$ep++;
			}

		}

		# Push any remaining groups onto the stack
		while ($startPos < $minHdrLen || $lsize > 0 || $esize > 0) {
			my $newStartPos = min($minHdrLen, $startPos + $progMaxRange);

			push @$lgroups, $lgroup;
			push @$esgroups, $esgroup;
			push @$edgroups, $edgroup;
			push @$shifts, $newStartPos - $startPos;

			$startPos = $newStartPos;
			$lgroup = [];
			$esgroup = [];
			$edgroup = [];
			$lsize = 0;
			$esize = 0;
		}

		$progLGroups{$hdrName} = $lgroups;
		$progESGroups{$hdrName} = $esgroups;
		$progEDGroups{$hdrName} = $edgroups;
		$progShifts{$hdrName} = $shifts;
	}

	# Work out the groups/shifts for each header type
	foreach my $hdrName (@hdrOrder) {
		my $hdr = $headers->{$hdrName};
		my $minHdrLen = $hdr->{'len'};
		my $numLens = 1;
		if ($minHdrLen == 0 && scalar(keys(%{$hdr->{'lenMap'}})) > 0) {
			$minHdrLen = min(values(%{$hdr->{'lenMap'}}));
			$numLens = scalar(keys(%{$hdr->{'lenMap'}}));
		}
		my $numNxtHdrs = scalar(keys(%{$hdr->{'nxtHdrMap'}}));
		$numNxtHdrs = 1 if $numNxtHdrs == 0;

		my $lookupCombs = $numLens * $numNxtHdrs;

		my $matchSeqs = [];
		my $maskSeqs = [];
		my $lens = [(0) x $lookupCombs];
		my $nxtHdrs = [(0) x $lookupCombs];
		for (my $i = 0; $i < $lookupCombs; $i++) {
			push @$matchSeqs, [(0) x $minHdrLen];
			push @$maskSeqs, [(0) x $minHdrLen];
		}

		my $nxtHdrMap = $hdr->{'nxtHdrMap'};
		my $lenMap = $hdr->{'lenMap'};


		my @hdrLens = ($hdr->{'len'});
		if ($hdrLens[0] == 0 && scalar(keys(%{$hdr->{'lenMap'}})) > 0) {
			@hdrLens = sort(cmpNum values(%{$hdr->{'lenMap'}}));
		}

		# Write the next header data
		my $pos = 0;
		if (scalar(keys(%$nxtHdrMap)) > 0) {
			my $nxtHdrBytes = $hdr->{'nxtHdrBytes'};
			foreach my $nxtHdrMatch (sort(keys(%$nxtHdrMap))) {
				my $nxtHdrName = $nxtHdrMap->{$nxtHdrMatch};
				my ($nhMatch, $nhMask) = extractMatchMaskBytes($nxtHdrMatch);
				for (my $i = 0; $i < $numLens; $i++) {
					for (my $j = 0; $j < scalar(@$nhMatch); $j++) {
						my $targ = $nxtHdrBytes->[$j];
						$matchSeqs->[$pos]->[$targ] |= $nhMatch->[$j];
						$maskSeqs->[$pos]->[$targ] |= $nhMask->[$j];
					}
					$nxtHdrs->[$pos] = $nxtHdrName;
					$pos++;
				}
			}
		}
		else {
			my $nxtHdrName = DONE;
			for (my $i = 0; $i < $numLens; $i++) {
				$nxtHdrs->[$pos] = $nxtHdrName;
				$pos++;
			}
		}

		# Write the length data
		$pos = 0;
		for (my $i = 0; $i < $numLens; $i++) {
			for (my $j = 0; $j < $numNxtHdrs; $j++) {
				my $len = $hdrLens[$i];
				if ($hdr->{'len'} == 0 && scalar(keys(%{$hdr->{'lenMap'}})) > 0) {
					my $lenBytes = $hdr->{'lenBytes'};
					my $lenMatch = findHashMatch($hdr->{'lenMap'}, $len);
					my ($lMatch, $lMask) = extractMatchMaskBytes($lenMatch);
					for (my $j = 0; $j < scalar(@$lMatch); $j++) {
						my $targ = $lenBytes->[$j];
						$matchSeqs->[$pos]->[$targ] |= $lMatch->[$j];
						$maskSeqs->[$pos]->[$targ] |= $lMask->[$j];
					}
				}
				$lens->[$pos] = $len;
				$pos++;
			}
		}

		$progMaskSeqs{$hdrName} = $maskSeqs;
		$progMatchSeqs{$hdrName} = $matchSeqs;
		$progLens{$hdrName} = $lens;
		$progNxtHdrs{$hdrName} = $nxtHdrs;
	}

	my $state = 0;
	my %states;
	my %entrySeen;

	# Add default DONE state
	addProgTableEntry(
		STATE => PROG_LOOKUP_DONE,
		SHIFT => $progMaxRange,
	);

	# Create state numbers for all header types
	foreach my $hdrName (@hdrOrder) {
		$states{$hdrName} = $state++;
	}
	$states{DONE} = PROG_LOOKUP_DONE;

	my $offsetIndex = 0;
	my %offsetMap = ();

	my %stateCount = ();

	# Finally, calculate the table entries
	foreach my $hdrName (@hdrOrder) {
		my $hdr = $headers->{$hdrName};
		my $minHdrLen = getMinHdrLen($hdr);

		my $lgroups = $progLGroups{$hdrName};
		my $esgroups = $progESGroups{$hdrName};
		my $edgroups = $progEDGroups{$hdrName};
		my $shifts = $progShifts{$hdrName};

		for (my $i = 0; $i < scalar(@{$progLens{$hdrName}}); $i++) {
			my $mask = $progMaskSeqs{$hdrName}->[$i];
			my $match = $progMatchSeqs{$hdrName}->[$i];
			my $len = $progLens{$hdrName}->[$i];
			my $nxtHdr = $progNxtHdrs{$hdrName}->[$i];

			my $hist = "$hdrName";
			my $currState = $states{$hist};
			my $startPos = 0;

			my $index = 0;
			while ($startPos < $len || $index < scalar(@$shifts)) {
				my $lgroup = $index < scalar(@$shifts) ? $lgroups->[$index] : [];
				my $esgroup = $index < scalar(@$shifts) ? $esgroups->[$index] : [];
				my $edgroup = $index < scalar(@$shifts) ? $edgroups->[$index] : [];
				my $shift;
				if ($index < scalar(@$shifts) - 1) {
					$shift = $shifts->[$index];
				} else {
					$shift = min($len - $startPos, $progMaxRange);
				}

				my @lMask = (0) x ($progLookupInputs * $progLookupWidth);
				my @lMatch = (0) x ($progLookupInputs * $progLookupWidth);

				# Need to populate the lookup mask/match if
				# values are being looked up
				if (scalar(@$lgroup) > 0) {
					for (my $lIndex = 0; $lIndex < scalar(@$lgroup); $lIndex++) {
						my $lPos = $lgroup->[$lIndex] + $startPos;
						for (my $k = $lPos; $k < $lPos + $progLookupWidth && $k < $minHdrLen; $k++) {
							my $dstIndex = $lIndex * $progLookupWidth + $k - $lPos;
							$lMask[$dstIndex] |= $mask->[$k];
							$lMatch[$dstIndex] |= $match->[$k];
						}
					}
				}

				my $lookupStr = maskMatchToLookupStr(\@lMask, \@lMatch);
				$hist .= ":$lookupStr";

				my $nextState;
				my $nxtHdrStr = $nxtHdr;
				$nxtHdrStr = $hist if ($startPos + $shift < $len);
				if (!defined($states{$nxtHdrStr})) {
					$states{$nxtHdrStr} = $state++;
				}
				$nextState = $states{$nxtHdrStr};

				if (!defined($entrySeen{$hist})) {
					$entrySeen{$hist} = 1;

					my $currOffsetIndex = 0;
					my $currOffsetAmt = 0;
					if ($progExtractHdrCnt{$hdrName} > 1) {
						if (!defined($offsetMap{$hdrName})) {
							$offsetMap{$hdrName} = ++$offsetIndex;
						}
						$currOffsetIndex = $offsetMap{$hdrName};
						$currOffsetAmt = $progExtractByteCnt{$hdrName};
					}

					addProgTableEntry(HIST => $hist,
						STATE => $currState,
						LOOKUP_STR => $lookupStr,
						LOOKUP_OFFSETS => $lgroup,
						EXTRACT_OFFSETS => $esgroup,
						EXTRACT_DSTS => $edgroup,
						SHIFT => $shift,
						NEXT_STATE => $nextState,
						OFFSET_INDEX => $currOffsetIndex,
						OFFSET_AMT => $currOffsetAmt);

					# Attempt to count unique states based
					# on location, match, and shift
					my $stateStr = "$hdrName:$startPos:$lookupStr:$shift";
					if (!defined($stateCount{$stateStr})) {
						$stateCount{$stateStr} = 0;
					}
					$stateCount{$stateStr}++;
				}

				$currState = $nextState;

				$startPos += $shift;
				$index++;
			}
		}

	}

	# Finally, walk through each of the states and grab the lookup fields
	# for the next type
	foreach my $state (sort(cmpNum keys(%progTableEntries))) {
		foreach my $entry (@{$progTableEntries{$state}}) {
			my $nextState = $entry->{NEXT_STATE};

			my $nextEntry = $progTableEntries{$nextState}->[0];
			my $nextEntryLGroup = $nextEntry->{LOOKUP_OFFSETS};

			$entry->{NEXT_LOOKUP_OFFSETS} = $nextEntryLGroup;
			$progNumEntries++;
		}
	}

	# Attempt to generate extra states for sequential types
	$progExtraEntries = 0;
	foreach my $hdrName (@hdrOrder) {
		#exit(1) if $hdrName eq 'icmpv6';
		#exit(1) if $hdrName eq 'ipv6';
		#exit(1) if $hdrName eq 'ipv4';
		#exit(1) if $hdrName eq 'eompls+ethernet2';
		my $currState = $states{$hdrName};

		# Work out how long the header is and what comes next
		print STDOUT "Considering: $hdrName ($currState)\n";
		my $pliRemaining = $progLookupInputs;
		my $lenRemaining = $progMaxRange;
		my $extRemaining = $progExtractWidth;
		my $currHdrName = $hdrName;
		my $elem = 0;
		my $depth = 0;
		my $startPos = 0;
		my @prevState;
		my @hdrSeqs;
		my @hdrSeq;
		my $hdrSeqStr = "";
		my %hdrSeqStrs;
		my %hdrCnts;
		#while ($pliRemaining > 0 && $lenRemaining > 0 && $currHdrName ne DONE && $elem < scalar(@{$progLens{$currHdrName}})) {
		while (1) {
			if ($elem >= scalar(@{$progLens{$currHdrName}})) {
				last;
			}
			#if ($pliRemaining > 0 && $lenRemaining > 0 && $currHdrName ne DONE && $elem < scalar(@{$progLens{$currHdrName}})) {

			#$hdrCnts{$currHdrName} = 0 if !defined($hdrCnts{$currHdrName});
			#print STDOUT "$currHdrName   SS: $hdrSeqStr   D: $depth   E: $elem   P: $pliRemaining   L: $lenRemaining   ER: $extRemaining   C: $hdrCnts{$currHdrName}   S: $startPos\n";
			my $currState = $states{$currHdrName};
			my $mask = $progMaskSeqs{$currHdrName}->[$elem];
			my $match = $progMatchSeqs{$currHdrName}->[$elem];
			my $len = $progLens{$currHdrName}->[$elem];
			my $nxtHdr = $progNxtHdrs{$currHdrName}->[$elem];
			my $maxHdrs = $progExtractHdrCnt{$currHdrName};
			$elem++;

			my $lgroups = $progLGroups{$currHdrName};
			my $esgroups = $progESGroups{$currHdrName};
			my $edgroups = $progEDGroups{$currHdrName};
			my $shifts = $progShifts{$currHdrName};

			my $pliDelta = 0;
			$pliDelta = scalar(@{$lgroups->[0]}) if defined($lgroups);
			my $extDelta = 0;
			$extDelta = scalar(@{$esgroups->[0]}) if defined($esgroups);

			my $prevPLIRemaining = $pliRemaining;
			my $prevLenRemaining = $lenRemaining;
			my $prevExtRemaining = $extRemaining;

			$pliRemaining -= $pliDelta;
			$lenRemaining -= $len;
			$extRemaining -= $extDelta;

			$hdrCnts{$currHdrName} = 0 if !defined($hdrCnts{$currHdrName});
			my $currHdrCnt = $hdrCnts{$currHdrName} + 1;

			#print STDOUT "--$currHdrName->$nxtHdr   SS: $hdrSeqStr   D: $depth   E: $elem   P: $pliRemaining   L: $lenRemaining   ER: $extRemaining   C: $hdrCnts{$currHdrName}   M: $maxHdrs\n";
			if ($startPos == 0 && $pliRemaining >= 0 && $currHdrCnt <= $maxHdrs && scalar(@$shifts) > 1) {
				push @hdrSeq, ($currHdrName, $len);
				$hdrSeqStr .= "$currHdrName:$len  ";

				if (scalar(@hdrSeq) > 2) {
					my @hdrSeqCopy = @hdrSeq;
					if (!defined($hdrSeqStrs{$hdrSeqStr})) {
						push @hdrSeqs, \@hdrSeqCopy;
						$hdrSeqStrs{$hdrSeqStr} = 1;
					}
				}

				pop @hdrSeq;
				pop @hdrSeq;
				$hdrSeqStr =~ s/\S+:\d+  $//;
			}

			if (scalar(@$shifts) == 1 && $lenRemaining >= 0 && $extRemaining >= 0 && $pliRemaining >= 0 && $currHdrCnt <= $maxHdrs) {
				push @hdrSeq, ($currHdrName, $len);
				$hdrSeqStr .= "$currHdrName:$len  ";

				if (scalar(@hdrSeq) > 2) {
					my @hdrSeqCopy = @hdrSeq;
					if (!defined($hdrSeqStrs{$hdrSeqStr})) {
						push @hdrSeqs, \@hdrSeqCopy;
						$hdrSeqStrs{$hdrSeqStr} = 1;
					}
				}

				if ($nxtHdr ne DONE && $lenRemaining >= 0) {
					$hdrCnts{$currHdrName}++;
					push @prevState, [$currHdrName, $elem, $prevPLIRemaining, $prevLenRemaining, $prevExtRemaining, $startPos];
					$currHdrName = $nxtHdr;
					$elem = 0;
					$startPos += $len;
					$depth++;
				} else {
					pop @hdrSeq;
					pop @hdrSeq;
					$hdrSeqStr =~ s/\S+:\d+  $//;
				}

			}

			if (scalar(@$shifts) > 1 || $lenRemaining < 0 || $pliRemaining < 0 || $extRemaining < 0 || $currHdrCnt > $maxHdrs || $nxtHdr eq DONE || $elem >= scalar(@{$progLens{$currHdrName}})) {
				$pliRemaining += $pliDelta;
				$lenRemaining += $len;
				$extRemaining += $extDelta;

				#if ($elem >= scalar(@{$progLens{$currHdrName}})){
				#	$hdrCnts{$currHdrName}--;
				#}

				while ($elem >= scalar(@{$progLens{$currHdrName}}) && $depth > 0) {

					my $prevState =  pop @prevState;
					($currHdrName, $elem, $pliRemaining, $lenRemaining, $extRemaining, $startPos) = @$prevState;
					$depth--;
					$hdrCnts{$currHdrName}--;

					pop @hdrSeq;
					pop @hdrSeq;
					$hdrSeqStr =~ s/\S+:\d+  $//;
				}
			}
		}

		# Remove subsets
		my %seqsByLen;
		for (my $i = 0; $i < scalar(@hdrSeqs); $i++) {
			my $hdrSeq = $hdrSeqs[$i];
			my $seqLen = scalar(@$hdrSeq) / 2;
			my $hdrSeqStr = join('::', @$hdrSeq);
			$seqsByLen{$seqLen} = {} if !defined($seqsByLen{$seqLen});
			$seqsByLen{$seqLen}{$hdrSeqStr} = $hdrSeq;
		}
		foreach my $seqLen (sort(keys(%seqsByLen))) {
			next if $seqLen == 1;
			next if !defined($seqsByLen{$seqLen - 1});

			foreach my $hdrSeqStr (keys(%{$seqsByLen{$seqLen}})) {
				my $shortStr = $hdrSeqStr;
				$shortStr =~ s/::[^:]+::[^:]+$//;

				if (defined($seqsByLen{$seqLen - 1}{$shortStr})) {
					delete $seqsByLen{$seqLen - 1}{$shortStr};
					#print STDOUT "  Deleting: $shortStr\n";
				}
			}
		}
		@hdrSeqs = ();
		foreach my $seqLen (sort(keys(%seqsByLen))) {
			foreach my $hdrSeqStr (sort(keys(%{$seqsByLen{$seqLen}}))) {
				my $hdrSeq = $seqsByLen{$seqLen}{$hdrSeqStr};
				push @hdrSeqs, $hdrSeq;
			}
		}

		# Identify the "best" sequence. Prefer:
		#   - most self-similar (ie. closest to the first header type)
		#   - shortest
		my $maxHdrs = 0;
		my $shortest = $progMaxRange + 1;
		my $chosen = -1;
		for (my $i = 0; $i < scalar(@hdrSeqs); $i++) {
			my $hdrSeq = $hdrSeqs[$i];
			my $seqLen = scalar(@$hdrSeq) / 2;
			my $hdrs = 0;
			for (my $j = 0; $j < scalar(@$hdrSeq); $j += 2) {
				$hdrs++ if ($hdrSeq->[$j] eq $hdrName);
			}
			if ($hdrs > $maxHdrs) {
				$maxHdrs = $hdrs;
				$shortest = $seqLen;
				$chosen = $i;
			} elsif ($hdrs == $maxHdrs && $seqLen < $shortest) {
				$shortest = $seqLen;
				$chosen = $i;
			}
			#print STDOUT "  Seq: " . join(' ', @$hdrSeq) . "\n";
		}
		if ($chosen >= 0) {
			print STDOUT "  Chosen: " . join(' ', @{$hdrSeqs[$chosen]}) . "\n";
			my @hdrSeq = @{$hdrSeqs[$chosen]};

			my @chosenSeqs;
			push @chosenSeqs, $hdrSeqs[$chosen];

			# Construct the combined group
			my @lgroup;
			my @esgroup;
			my @edgroup;
			my $shift = 0;
			for (my $i = 0; $i < scalar(@hdrSeq); $i += 2) {
				my $hName = $hdrSeq[$i];
				my $hLen = $hdrSeq[$i + 1];

				my @hLGroup = @{$progLGroups{$hName}->[0]};
				my @hESGroup = @{$progESGroups{$hName}->[0]};
				my @hEDGroup = @{$progEDGroups{$hName}->[0]};

				@hLGroup = map { $_ + $shift; } @hLGroup;
				@hESGroup = map { $_ + $shift; } @hESGroup;

				push @lgroup, @hLGroup;
				push @esgroup, @hESGroup;
				push @edgroup, @hEDGroup;
				$shift += $hLen;
			}

			#print STDOUT "    LGroup: [@lgroup]   ESGroup: [@esgroup]   EDGroup: [@edgroup]   Shift: $shift\n";

			# Identify if any of the other groups match
			my %otherSeqStrs;
			for (my $i = 0; $i < scalar(@hdrSeqs); $i++) {
				next if $i == $chosen;

				my $hdrSeq = $hdrSeqs[$i];
				my @otherSeq;
				my $otherSeqStr = "";
				my $base = 0;
				my $lookupPos = 0;
				my $done = 0;
				for (my $j = 0; $j < scalar(@$hdrSeq) && !$done; $j += 2) {
					my $hName = $hdrSeq->[$j];
					my $hLen = $hdrSeq->[$j + 1];

					my @hLGroup = @{$progLGroups{$hName}->[0]};
					@hLGroup = map { $_ + $base; } @hLGroup;
					$base += $hLen;

					$done |= scalar(@hLGroup) + $lookupPos > scalar(@lgroup);

					for (my $k = 0; $k < scalar(@hLGroup) && !$done; $k++) {
						$done = $hLGroup[$k] != $lgroup[$lookupPos++];
					}

					if (!$done) {
						push @otherSeq, ($hName, $hLen);
						$otherSeqStr .= "$hName:$hLen  ";
					}
				}

				if (scalar(@otherSeq) > 2) {
					$otherSeqStrs{$otherSeqStr} = [] if !defined($otherSeqStrs{$otherSeqStr});
					push @{$otherSeqStrs{$otherSeqStr}}, $i;
					push @chosenSeqs, \@otherSeq;
				}
			}

			for (my $i = 1; $i < scalar(@chosenSeqs); $i++) {
				my $hdrSeq = $chosenSeqs[$i];
				#print STDOUT "  Other seq: " . join(' ', @$hdrSeq) . "\n";
			}

			my %extraStateStrs;
			for (my $i = 0; $i < scalar(@chosenSeqs); $i++) {
				my $hdrSeq = $chosenSeqs[$i];
				my $seqStr = "";
				my @seq;
				for (my $j = 0; $j < scalar(@$hdrSeq); $j += 2) {
					my $hName = $hdrSeq->[$j];
					my $hLen = $hdrSeq->[$j + 1];

					push @seq, $hName, $hLen;
					$seqStr .= "$hName";

					if ($j >= 2) {
						if (!defined($extraStateStrs{$seqStr})) {
							$extraStateStrs{$seqStr} = 1;
							$progExtraEntries++;
							my $hdr = $headers->{$hName};
							foreach my $nxtState (values(%{$hdr->{'nxtHdrMap'}})) {
								my $nxtSeqStr = $seqStr . ":$hLen  $nxtState";
								if (!defined($extraStateStrs{$nxtSeqStr})) {
									$extraStateStrs{$nxtSeqStr} = 1;
									$progExtraEntries++;
								}
							}
						}
					}
					$seqStr .= ":$hLen  ";
				}
			}
			foreach my $extraState (keys(%extraStateStrs)) {
				print STDOUT "  Extra state: $extraState\n";
			}
		}
	}


	# Count the number of unique entries
	foreach my $stateStr (sort(keys(%stateCount))) {
		my $cnt = $stateCount{$stateStr};
		$progUniqueEntries++;
	}

	# Get the lookup fields for the first state
	$progFirstLookupWords = $progTableEntries{0}->[0]->{LOOKUP_OFFSETS};

	return ($progNumEntries, $progUniqueEntries, $progExtraEntries, \%progTableEntries, $progFirstLookupWords);
}

sub getProgLookupTable {
	return ($progNumEntries, $progUniqueEntries, $progExtraEntries, \%progTableEntries, $progFirstLookupWords);
}

# Convert a mask and match array into a programmable parser table lookup string
sub maskMatchToLookupStr {
	my ($mask, $match) = @_;

	my $ret = "";
	for (my $i = 0; $i < $progLookupInputs; $i++) {
		$ret .= ', ' if $i > 0;

		$ret .= (8 * $progLookupWidth) . "'b";

		for (my $j = 0; $j < $progLookupWidth; $j++) {
			my $maskB = $mask->[$i * $progLookupWidth + $j];
			my $matchB = $match->[$i * $progLookupWidth + $j];

			for (my $bit = 0; $bit < 8; $bit++) {
				$ret .= '_' if ($bit % 4 == 0);
				if (($maskB & (2 ** (7 - $bit))) == 0) {
					$ret .= '?';
				}
				else {
					$ret .= chr(ord('0') + (($matchB >> (7 - $bit)) & 1));
				}
			}
		}
	}

	return $ret;
}

# Add a programmable table entry
sub addProgTableEntry {
	my %args = @_;

	# Verify we have a state
	confess "Must specify a state" unless defined($args{STATE});

	# Set defaults
	$args{LOOKUP_OFFSETS} = [] unless defined($args{LOOKUP_OFFSETS});
	$args{EXTRACT_OFFSETS} = [] unless defined($args{EXTRACT_OFFSETS});
	$args{EXTRACT_DSTS} = [] unless defined($args{EXTRACT_DSTS});
	$args{NEXT_LOOKUP_OFFSETS} = [] unless defined($args{NEXT_LOOKUP_OFFSETS});
	$args{NEXT_STATE} = PROG_LOOKUP_DONE unless defined($args{NEXT_STATE});
	if (!defined($args{LOOKUP_STR})) {
		my $maskLen = $progLookupWidth * $progLookupInputs;
		$args{LOOKUP_STR} = maskMatchToLookupStr([(0) x $maskLen],
			[(0) x $maskLen]);
	}
	foreach my $progField (@allProgFields) {
		$args{$progField} = 0 unless defined($args{$progField});
	}

	# Add the entry
	$progTableEntries{$args{STATE}} = []
		unless defined($progTableEntries{$args{STATE}});

	push @{$progTableEntries{$args{STATE}}}, \%args;
}

sub setProgExtractWidth {
	$progExtractWidth = shift;
}

sub setProgLookupInputs {
	$progLookupInputs = shift;
}

sub setProgStateWidth {
	$progStateWidth = shift;
}

sub setProgLookupWidth {
	$progLookupWidth = shift;
}

sub setProgMaxRange {
	$progMaxRange = shift;
}

sub setProgFirstLookupWordsNew {
	$progFirstLookupWordsNew = shift;
}

sub getProgFirstLookupWordsNew {
	return $progFirstLookupWordsNew;
}

sub expandFilename {
	my $filename = shift;
	$filename =~ s{ ^ ~ ( [^/]* ) }
		      { $1
			      ? (getpwnam($1))[7]
			      : ( $ENV{HOME} || $ENV{LOGDIR}
				    || (getpwuid($<))[7]
				)
	  }ex;

	return $filename;
}

# Convert a mask and match array into a programmable parser table lookup string
sub maskMatchWithStateToLookupStr {
	my ($mask, $match) = @_;

	my $ret = "";
	my $stateBytes = POSIX::ceil($progStateWidth / 8.0);
	for (my $i = 0; $i < $stateBytes; $i++) {
		my $maskB = $mask->[$i];
		my $matchB = $match->[$i];
		for (my $bit = 0; $bit < 8; $bit++) {
			$ret .= '_' if ($bit % 4 == 0);
			if (($maskB & (2 ** (7 - $bit))) == 0) {
				$ret .= '?';
			}
			else {
				$ret .= chr(ord('0') + (($matchB >> (7 - $bit)) & 1));
			}
		}
	}
	$ret = "${progStateWidth}'b" . substr($ret, $stateBytes * 10 - $progStateWidth - $progStateWidth / 4);

	for (my $i = 0; $i < $progLookupInputs; $i++) {
		$ret .= ', ';

		$ret .= (8 * $progLookupWidth) . "'b";

		for (my $j = 0; $j < $progLookupWidth; $j++) {
			my $maskB = $mask->[$i * $progLookupWidth + $j + $stateBytes];
			my $matchB = $match->[$i * $progLookupWidth + $j + $stateBytes];

			for (my $bit = 0; $bit < 8; $bit++) {
				$ret .= '_' if ($bit % 4 == 0);
				if (($maskB & (2 ** (7 - $bit))) == 0) {
					$ret .= '?';
				}
				else {
					$ret .= chr(ord('0') + (($matchB >> (7 - $bit)) & 1));
				}
			}
		}
	}

	return $ret;
}

sub parseTCAMEntry {
	my $line = shift;
	
	my $entry = {};
	if ($line =~ /Match: \(\[([^\]]*)\], \[([^\]]*)\]\)/) {
		my @mask = split(', ', $1);
		@mask = map(hex, @mask);

		my @match = split(', ', $2);
		@match = map(hex, @match);

		my $lookupStr = maskMatchWithStateToLookupStr(\@mask, \@match);

		$entry->{'MASK'} = \@mask;
		$entry->{'MATCH'} = \@match;
		$entry->{'LOOKUP_STR'} = $lookupStr;
	}
	if ($line =~ /Next-State:\s*(\d+)\/(\d+)/) {
		$entry->{'NXT_STATE_VAL'} = int($1);
		$entry->{'NXT_STATE_MASK'} = int($2);
	}
	if ($line =~ /Adv:\s*(\d+)/) {
		$entry->{'SHIFT'} = $1;
	}
	if ($line =~ /Next-Lookup:\s*\[([^\]]*)\]/) {
		my @nxtLookups = map(int, split(', ', $1));

		$entry->{'LOOKUP_OFFSETS'} = \@nxtLookups;
	}
	if ($line =~ /Hdr-Starts:\s*\[([^\]]*)\]/) {
		my $hdrStarts = $1;
		$hdrStarts =~ s/\s|\(|\)//g;
		my @hdrStarts = split(',', $hdrStarts);

		$entry->{'HDR_STARTS'} = \@hdrStarts;
	}
	if ($line =~ /Extract:\s*\[([^\]]*)\]/) {
		my $extract = $1;
		$extract =~ s/\s|\(|\)//g;
		my @extract = split(',', $extract);

		$entry->{'EXTRACT'} = \@extract;
	}

	if ($line =~ /# Match:\s*(\[.*\])\s+Nxt-State:\s(\d+)/) {
		$entry->{'MATCH_CMT'} = $1;
		$entry->{'NXT_STATE_CMT'} = $2;
	}

	return $entry;
}

sub parseFirstLookupEntry {
	my $line = shift;
	
	my $entry = {};
	if ($line =~ /First-Lookup:\s*\[([^\]]*)\]/) {
		my @nxtLookups = map(int, split(', ', $1));

		$entry->{'LOOKUP_OFFSETS'} = \@nxtLookups;
	}

	return $entry;
}

sub setTCAMTableFile {
	$tcamTableFile = shift;
}

sub getTCAMTableFile {
	return $tcamTableFile;
}

1;
