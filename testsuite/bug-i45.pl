#! /usr/bin/perl

# bug i45 - cmsearch incorrect model boundaries in alidisplay with local ends at terminii
#
# EPN, Tue Feb 21 10:27:18 2017
#

$usage = "perl bug-i45.pl <cmsearch> <path to bug-i45.cm> <path to bug-i45.fa>\n";
if ($#ARGV != 1) { die "Wrong argument number.\n$usage"; }

$cmsearch = shift;
$cmfile   = shift;
$seqfile  = shift;
$ok       = 1;

# Model is RF00177 from Rfam 12.1
# Sequence is gi|1029367024|emb|FKPB01000027.1|
# This is the original CM/sequence pair that exposed the bug to Azat Badretdin.
#
# Original bug (latest version present v1.1.2) causes cm_alidisplay to have a
# model boundary on the 3' end (2367) that exceeds the model's consensus length (1533).
# 
# This test script checks that the 5' model boundary in the alignment 
# is 1 and that th e3' model boundary in the alignment is 1533, 
# and fails if they are not.
#
if ($ok) { 
  $output = `$cmsearch --cpu 0 --onepass --5trunc --rfam --notextw $cmfile $seqfile`
  if ($? != 0) { $ok = 0; }
  # find the relevant line of the output, and make sure it has 1 as the first
  # model position (2nd token) and 1533 as the final one (4th token)
  if($output !~ /\n\s+SSU\_rRNA\_bacteria\s+1\s+\S+\s+1533\n/) { 
    $ok = 0;
  }
}

if ($ok) { print "ok\n";     exit 0; }
else     { print "FAILED\n"; exit 1; }
