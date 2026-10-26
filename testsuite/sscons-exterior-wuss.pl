#! /usr/bin/perl

# Consensus structure lines written by Infernal (cmbuild -O, cmalign,
# cmemit -a) must use full WUSS notation. Single-stranded consensus
# positions that are not enclosed by any base pair (external) should
# be annotated ':'. Infernal 1.1.5 and earlier annotated external
# positions that lie between two top-level stems as '-' (bulge/interior)
# or ',' (multiloop) instead, because CreateCMConsensus() decided
# external-vs-internal from a face count that included the stems of
# the sibling subtree of a top-level BIF node.
#
# Two small alignments are used: one with two top-level stems, and
# one with three top-level stems, an interior loop and a multiloop
# (so that '-' and ',' are still required inside closed loops). The
# expected strings are full WUSS, matching Easel's esl_wuss_full()
# conversion of the input structure.
#
# EPN, Mon Oct 26 2026 (w/Claude)
#
$usage = "perl sscons-exterior-wuss.pl <cmbuild> <cmalign> <cmemit>\n";
if ($#ARGV != 2) { die "Wrong argument number.\n$usage"; }

$cmbuild = shift;
$cmalign = shift;
$cmemit  = shift;
$ok      = 1;

@ssA = ("::<<<<____>>>>:::<<<<____>>>>::",
        "::<<<<___>>>>::<<<--<<<___>>>>>>::(((,,<<<___>>>,,<<<___>>>,,)))::");

@tmpfileA = ();
for ($a = 0; $a < scalar(@ssA) && $ok; $a++) {
  $ss  = $ssA[$a];
  $seq = make_seq($ss);
  $sto = "sscons-ext.$a.sto";
  $cm  = "sscons-ext.$a.cm";
  $bld = "sscons-ext.$a.O.stk";
  $aln = "sscons-ext.$a.aln.stk";
  $emt = "sscons-ext.$a.emit.stk";
  push(@tmpfileA, $sto, $cm, $bld, $aln, $emt);

  open(OUT, ">$sto") || die;
  print OUT "# STOCKHOLM 1.0\n\n";
  print OUT "seq1          $seq\n";
  print OUT "seq2          $seq\n";
  print OUT "seq3          $seq\n";
  print OUT "#=GC SS_cons  $ss\n";
  print OUT "//\n";
  close(OUT);

  if ($ok) { system("$cmbuild -F -O $bld $cm $sto > /dev/null 2>&1");           if ($? != 0) { $ok = 0; } }
  if ($ok) { system("$cmalign -o $aln $cm $sto > /dev/null 2>&1");               if ($? != 0) { $ok = 0; } }
  if ($ok) { system("$cmemit -a -N 3 --seed 7 -o $emt $cm > /dev/null 2>&1");    if ($? != 0) { $ok = 0; } }

  foreach $file ($bld, $aln, $emt) {
    if ($ok) {
      $got = consensus_ss($file);
      if ($got ne $ss) {
        printf("SS_cons mismatch in $file\n  expected: $ss\n  got:      $got\n");
        $ok = 0;
      }
    }
  }
}

foreach $tmpfile (@tmpfileA) {
  unlink $tmpfile if -e $tmpfile;
}

if ($ok) { print "ok\n";     exit 0; }
else     { print "FAILED\n"; exit 1; }


# make_seq(): a sequence that matches <ss>: G-C base pairs, A elsewhere.
sub make_seq {
  my ($ss) = @_;
  my @ssA  = split("", $ss);
  my @seqA = ("A") x scalar(@ssA);
  my @stk  = ();
  for (my $i = 0; $i < scalar(@ssA); $i++) {
    if    ($ssA[$i] =~ m/[\<\(\[\{]/) { push(@stk, $i); }
    elsif ($ssA[$i] =~ m/[\>\)\]\}]/) { my $j = pop(@stk); $seqA[$j] = "G"; $seqA[$i] = "C"; }
  }
  return join("", @seqA);
}

# consensus_ss(): the SS_cons of a (single block) Stockholm file, with
# insert ('.') and local-end ('~') padding columns removed.
sub consensus_ss {
  my ($file) = @_;
  my $ss = "";
  open(IN, $file) || die "unable to open $file";
  while ($line = <IN>) {
    if ($line =~ m/^\#=GC\s+SS_cons\s+(\S+)/) { $ss .= $1; }
  }
  close(IN);
  $ss =~ s/[\.\~]//g;
  return $ss;
}
