#! /usr/bin/perl

# bug i64 - HMM-banded Outside writes local-end (EL) cells at d < 0
#
# In cm_OutsideAlignHB()/cm_TrOutsideAlignHB() each local-end-capable
# state v adds into the EL deck at d = (v's d) - sd. For an emitting
# state whose band reaches d < sd that is d < 0, i.e. the last cell of
# the previous EL row (in bounds, so memory checkers are silent). The
# stray value changes posteriors and can change the optimal-accuracy
# alignment of default (local, truncation-aware) cmalign; bit scores
# are unaffected.
#
# Test model: bug-i64.sto is the Rfam 15.10 RF02194 (HPnc0260)
# seed alignment with all #=GF/#=GS/#=GR annotation removed, built with
# 'cmbuild -F' as Rfam does. Target: seed sequence
# AE017125.1/1703851-1703967. With the bug, default cmalign puts one
# 39-consensus-position local end in the second hairpin; fixed, it puts
# two local ends (20 positions in the first hairpin, 19 in the second).
# We check that the #=GC SS_cons line has two runs of '~' (local ends),
# not one.
#
# EPN, Fri Oct  9 2026 (w/Claude)
#

$usage = "perl bug-i64.pl <cmbuild> <cmalign> <testsuite/bug-i64.sto>\n";
if ($#ARGV != 2) { die "Wrong argument number.\n$usage"; }

$cmbuild = shift;
$cmalign = shift;
$alifile = shift;
$ok      = 1;

# Make our test sequence to align, i64.1
#
open (OUT, ">i64.1") || die;
print OUT <<END;
>AE017125.1/1703851-1703967
AAUUUUUAACCUUGCUGAUGUGCUCAUUGAUAUAAGUGUGGUGUUGAUGAUUUAUCAAAU
GUUACAAGAAAAAAGAAAUGGUGGCCACGACUGGACUUGAACCAGCGACCACUACCA
END
close OUT;

# build the model
if ($ok) {
  $output = `$cmbuild -F i64.cm $alifile`;
  if ($? != 0) { $ok = 0; }
  # make sure output includes a 'CPU' string indicating success
  if($output !~ /\n\# CPU/) {
    $ok = 0;
  }
}

# align with default options (local, truncation-aware), one line per seq
if ($ok) {
  $output = `$cmalign --outformat pfam -o i64.2 i64.cm i64.1`;
  if ($? != 0) { $ok = 0; }
  # make sure output includes a 'CPU' string indicating success
  if($output !~ /\n\# CPU/) {
    $ok = 0;
  }
}

# count the local ends ('~' runs) in SS_cons
if ($ok) {
  $nel = -1;
  if(open(IN, "i64.2")) {
    while($line = <IN>) {
      if($line =~ /^\#=GC SS_cons\s+(\S+)/) {
        $sscons = $1;
        $nel = 0;
        while($sscons =~ /~+/g) { $nel++; }
      }
    }
    close(IN);
  }
  if($nel != 2) { $ok = 0; }
}

foreach $tmpfile ("i64.1", "i64.2", "i64.cm") {
  unlink $tmpfile if -e $tmpfile;
}
if ($ok) { print "ok\n";     exit 0; }
else     { print "FAILED\n"; exit 1; }
