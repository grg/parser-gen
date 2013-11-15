#!/usr/bin/perl

use strict ;
use warnings ;
use diagnostics -verbose ;

my $log_file = shift or die "No log file supplied as arguments" ;
-e $log_file or die "No such logfile exists\n" ;

my $status = `cat $log_file | grep -v SCRIPT-Error` ;

print "++++++++++++++++++++++++++++++++++\n" ;
print "Checking $log_file in checkRun.pl\n" ;
print "++++++++++++++++++++++++++++++++++\n" ;

foreach $status (split('\n', $status)) {
	$status =~ m/Fatal:/    and die "Fatal error in run\n" ;
	$status =~ m/ERROR/     and die "ERROR in run\n" ;
	$status =~ m/Error:/    and die "Error in run\n" ;
	$status =~ m/Abort at / and die "Abort in run\n" ;
	$status =~ m/^Error-\[.*\]/ and die "Error in run\n" ;
}
