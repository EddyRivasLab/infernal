#! /usr/bin/perl

# bug i52 - cmalign --nonbanded SIGSEGVs on truncated alignment of some models
#
# The d==1 case of the MP emission block is hand-unrolled out of the
# 'for (d = 2; d <= j; d++)' loop that is the only thing enforcing d <= j, so it
# also runs at j==0. There i = j-d+1 = 0, and lmesc[v][dsq[0]] / rmesc[v][dsq[0]]
# index the digitized-sequence sentinel byte (255) instead of a residue, reading
# 237 elements past the state's row of a flat Kp*nstates allocation.
#
# The process only crashes when the float at that address happens to be a NaN
# bit pattern that then reaches FLogsum, which is why it is model-dependent.
#
# *** NOTE ON WHAT THIS TEST CAN AND CANNOT DETECT ***
# Because the crash depends on what sits past the end of the allocation, it is
# HEAP-LAYOUT DEPENDENT. Both models below segfault on every released version
# back to 1.1.2 and on develop before the fix, but a future allocator, platform
# or compiler could place a non-NaN float there, in which case these runs would
# complete on UNFIXED code and this test would silently stop discriminating.
# If you are auditing this test, the check that it still has teeth is: build a
# binary without the j==0 guards in cm_dpalign_trunc.c and confirm both models
# still segfault. Two independent models are used to make that less likely to
# happen to both at once.
#
# The .sto files are the alignments the .cm files were built from, committed as
# provenance only -- this test does not read them. Regenerate with:
#   cmbuild -F bug-i52.cm bug-i52.sto
# The committed .cm files are what the test runs, deliberately: rebuilding here
# would make the test vulnerable to unrelated cmbuild changes.
#
# EPN, Mon Sep 14 2026
#

$usage = "perl bug-i52.pl <cmalign> <path to bug-i52.cm> <path to bug-i52.fa> <path to bug-i52b.cm> <path to bug-i52b.fa>\n";
if ($#ARGV != 4) { die "Wrong argument number.\n$usage"; }

$cmalign  = shift;
$cmfile1  = shift;
$seqfile1 = shift;
$cmfile2  = shift;
$seqfile2 = shift;
$ok       = 1;

# Both objectives, both models. Before the fix every one of these segfaults.
foreach $pair ("$cmfile1 $seqfile1", "$cmfile2 $seqfile2") {
  foreach $opt ("--optacc", "--cyk") {
    if ($ok) {
      $output = `$cmalign --cpu 0 --nonbanded $opt -o i52.sto $pair > /dev/null 2>&1`;
      if ($? != 0) { $ok = 0; }
      # a signal death leaves no output file, or an empty one
      if ($ok) {
        if   (! -e "i52.sto") { $ok = 0; }
        elsif(! -s "i52.sto") { $ok = 0; }
      }
    }
  }
}

foreach $tmpfile ("i52.sto") {
  unlink $tmpfile if -e $tmpfile;
}

if ($ok) { print "ok\n";     exit 0; }
else     { print "FAILED\n"; exit 1; }
