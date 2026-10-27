#! /usr/bin/perl

# bug: cmbuild entropy weighting fails with "xl,xr do not bracket a root"
#
# cm_EntropyWeight() finds the effective sequence number (Neff) that
# gives the target relative entropy by bisection. After checking
# that f(min_Neff) < 0 < f(max_Neff), where f is relative entropy
# minus the target, it bisected over [0, max_Neff], not [min_Neff,
# max_Neff]. Relative entropy is not monotonic in Neff near 0: with
# the default priors, the base pair relative entropy of the prior
# alone (Neff = 0) can exceed that at small Neff. For a model with
# many base pairs and a low enough target, f(0) >= 0 and the
# rootfinder failed (abort) because [0, max_Neff] does not bracket a
# root. Seen with --ere 0.38 on some Rfam seed alignments, and with
# default settings on an almost entirely base-paired alignment.
#
# bug-eweight-bracket.sto is a synthetic alignment (random base
# pairs, 8 sequences, 60 base pairs and a 4 nt loop). Its relative
# entropy is about 0.572 bits at Neff 0 and 0.554 at Neff 0.1, so
# '--ere 0.565' failed. We check that cmbuild succeeds and that the
# effective sequence number is between 0.1 (the default --eminseq)
# and the number of sequences.
#
# EPN, Tue Oct 27 2026 (w/Claude)
#

$usage = "perl bug-eweight-bracket.pl <cmbuild> <testsuite/bug-eweight-bracket.sto>\n";
if ($#ARGV != 1) { die "Wrong argument number.\n$usage"; }

$cmbuild = shift;
$alifile = shift;
$ok      = 1;
$cmfile  = "bug-eweight-bracket.cm";

$output = `$cmbuild --ere 0.565 --p7eent -F $cmfile $alifile 2>&1`;
if ($? != 0)                  { $ok = 0; }
if ($output !~ /\n\# CPU/)    { $ok = 0; }

# check the effective sequence number
$effn = -1.;
if ($ok) {
  if (open(CM, $cmfile)) {
    while ($line = <CM>) {
      if ($line =~ /^EFFN\s+(\S+)/) { $effn = $1; last; }
    }
    close(CM);
  }
  if ($effn < 0.1 || $effn > 8.) { $ok = 0; }
}

unlink $cmfile if -e $cmfile;
if ($ok) { print "ok\n";     exit 0; }
else     { print "FAILED\n"; exit 1; }
