#! /usr/bin/perl

# bug: a CM longer than p7_MAXM aborts every program with a HMMER
# fatal exception
#
# Every CM carries profile HMMs as long as its consensus, and HMMER
# limits profile HMM length to p7_MAXM, a build-time setting.
# Reading a CM with a longer consensus aborted (SIGABRT) with
# "Fatal exception (source file p7_hmm.c ...): model length M=...
# out of range", naming neither the model nor the limit. The CM
# readers now check CLEN against p7_MAXM first and fail with a
# normal error.
#
# We set the CLEN of the second model in 4.c.cm to 2,000,000,000
# (beyond any p7_MAXM) and check that cmstat reports the first model,
# then exits with status 1 and an error that names the second model
# and p7_MAXM.
#
# EPN, Wed Oct 28 2026 (w/Claude)
#

$usage = "perl bug-clen-maxm.pl <cmstat> <testsuite/4.c.cm>\n";
if ($#ARGV != 1) { die "Wrong argument number.\n$usage"; }

$cmstat = shift;
$cmfile = shift;
$ok     = 1;
$tmpcm  = "bug-clen-maxm.cm";

# copy, changing the second model's CLEN; remember its name
open(IN,  $cmfile)   || die "FAILED: couldn't open $cmfile";
open(OUT, ">$tmpcm") || die "FAILED: couldn't open $tmpcm for writing";
$nmodel = 0;
$name2  = "";
while ($line = <IN>) {
  if ($line =~ /^INFERNAL1/) { $nmodel++; }
  if ($nmodel == 2 && $line =~ /^NAME\s+(\S+)/) { $name2 = $1; }
  if ($nmodel == 2 && $line =~ /^CLEN\s+/)      { $line = "CLEN     2000000000\n"; }
  print OUT $line;
}
close(IN);
close(OUT);
if ($name2 eq "") { $ok = 0; }

$output = `$cmstat $tmpcm 2>&1`;
$status = $?;
if ($status != (1 << 8))                                                     { $ok = 0; } # exit status 1, not a signal
if ($output !~ /model \Q$name2\E: consensus length 2000000000 exceeds p7_MAXM/) { $ok = 0; }
if ($output !~ /\n\s*1\s+\S+\s+RF00005/)                                     { $ok = 0; } # first model still reported

unlink $tmpcm if -e $tmpcm;
if ($ok) { print "ok\n";     exit 0; }
else     { print "FAILED\n"; print $output; exit 1; }
