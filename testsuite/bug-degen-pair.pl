#! /usr/bin/perl

# bug: ParsetreeScore() read out of bounds when computing the
# structure score of a base pair that includes a degenerate residue
#
# For a base pair with a degenerate residue (e.g. R or N),
# ParsetreeScore() added the correct (averaged) pair score to the
# total score but indexed the pair score table with the degenerate
# symbol for the structure score, reading past the end of the
# table. The "STRUCTURE SCORE" lines written by the --tfile options
# of cmbuild, cmalign and cmemit were wrong for such sequences.
#
# bug-degen-pair.sto is a synthetic hairpin alignment. Rows seqA,
# seqG and seqR are identical except at one paired column, where
# they have A, G and R (A or G), each paired with U. With the default
# uniform background, the degenerate R-U pair scores the mean of the
# A-U and G-U pair scores, and the marginal scores of R are the mean
# of those of A and G. So both the total score and the structure
# score of seqR must be the mean of those of seqA and seqG. We build
# with 'cmbuild --tfile' and check both, to the 0.01 bit precision
# of the --tfile output.
#
# EPN, Mon Oct 19 2026 (w/Claude)
#

$usage = "perl bug-degen-pair.pl <cmbuild> <testsuite/bug-degen-pair.sto>\n";
if ($#ARGV != 1) { die "Wrong argument number.\n$usage"; }

$cmbuild = shift;
$alifile = shift;
$ok      = 1;
$tol     = 0.015; # scores are printed with 2 decimals

if ($ok) {
  $output = `$cmbuild -F --tfile degen-pair.tfile degen-pair.cm $alifile`;
  if ($? != 0) { $ok = 0; }
  # make sure output includes a 'CPU' string indicating success
  if($output !~ /\n\# CPU/) {
    $ok = 0;
  }
}

# parse SCORE and STRUCTURE SCORE for each sequence
if ($ok) {
  %scH     = ();
  %structH = ();
  $sqname  = "";
  if (open(TF, "degen-pair.tfile")) {
    while($line = <TF>) {
      if   ($line =~ /^>(\S+)/)                                 { $sqname = $1; }
      elsif($line =~ /^\s+STRUCTURE SCORE:\s+(\S+) bits/)       { $structH{$sqname} = $1; }
      elsif($line =~ /^\s+SCORE:\s+(\S+) bits/)                 { $scH{$sqname}     = $1; }
    }
    close(TF);
  }
  else { $ok = 0; }
}

if ($ok) {
  foreach $sqname ("seqA", "seqG", "seqR") {
    if((! exists $scH{$sqname}) || (! exists $structH{$sqname})) { $ok = 0; }
  }
}

if ($ok) {
  if(abs($scH{"seqR"}     - (($scH{"seqA"}     + $scH{"seqG"})     / 2.)) > $tol) { $ok = 0; }
  if(abs($structH{"seqR"} - (($structH{"seqA"} + $structH{"seqG"}) / 2.)) > $tol) { $ok = 0; }
}

foreach $tmpfile ("degen-pair.tfile", "degen-pair.cm") {
  unlink $tmpfile if -e $tmpfile;
}
if ($ok) { print "ok\n";     exit 0; }
else     { print "FAILED\n"; exit 1; }
