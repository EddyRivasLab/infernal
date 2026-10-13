#! /usr/bin/perl

# bug i66 - cmbuild -O corrupts truncated fragment rows of zero-basepair models
#
# For models with zero basepairs, cmbuild runs cm_parsetree_Doctor()
# on each training parsetree. For a 5'-truncated fragment row whose
# parsetree starts with a run of marginal (non-left-emitting) nodes,
# e.g. because the row has residues in the 3' flanking insert, Doctor
# forced those nodes to emit, on missing-data ('~') columns. The CM
# was unaffected, but the -O alignment writer read outside the
# sequence for them and wrote garbage bytes or dropped the row's
# residues (on Linux usually silently: a well-formed row of all '~').
# The same Doctor check made cmbuild fail outright on a fragment
# truncated at both ends whose first residue is in an insert column.
#
# bug-i66.sto is a synthetic (random sequence) file with two zero-bp
# alignments. In i66a, frag1-frag3 are 5'-truncated fragments with
# 3' flanking insert residues (frag3 also has a deletion). In i66b,
# frag4 is truncated at both ends and starts in an insert column. We
# build with 'cmbuild --hand --noss -F -O' and check that cmbuild
# succeeds, that the -O file has only printable ASCII, and that every
# sequence's residues in the -O file match the input.
#
# EPN, Tue Oct 13 2026 (w/Claude)
#

$usage = "perl bug-i66.pl <cmbuild> <testsuite/bug-i66.sto>\n";
if ($#ARGV != 1) { die "Wrong argument number.\n$usage"; }

$cmbuild = shift;
$alifile = shift;
$ok      = 1;

# build the models, saving the -O alignments
if ($ok) {
  $output = `$cmbuild --hand --noss -F -O i66.1 i66.cm $alifile`;
  if ($? != 0) { $ok = 0; }
  # make sure output includes a 'CPU' string indicating success
  if($output !~ /\n\# CPU/) {
    $ok = 0;
  }
}

# compare residues of each sequence in each alignment, input vs -O output
if ($ok) {
  @in  = read_stockholm_residues($alifile);
  @out = read_stockholm_residues("i66.1");
  if(scalar(@in) != 2 || scalar(@out) != 2) { $ok = 0; }
  for($a = 0; $ok && $a < scalar(@in); $a++) {
    %inH  = %{$in[$a]};
    %outH = %{$out[$a]};
    if(scalar(keys %inH) != scalar(keys %outH)) { $ok = 0; }
    foreach $sqname (keys %inH) {
      if((! exists $outH{$sqname}) || ($inH{$sqname} ne $outH{$sqname})) { $ok = 0; }
    }
  }
}

# -O output must contain only printable ASCII (plus newlines)
if ($ok) {
  if(open(IN, "i66.1")) {
    while($line = <IN>) {
      if($line =~ /[^\x20-\x7e\n]/) { $ok = 0; }
    }
    close(IN);
  }
  else { $ok = 0; }
}

foreach $tmpfile ("i66.1", "i66.cm") {
  unlink $tmpfile if -e $tmpfile;
}
if ($ok) { print "ok\n";     exit 0; }
else     { print "FAILED\n"; exit 1; }

# read_stockholm_residues(<file>): return a list with one hash
# reference per alignment, mapping each sequence name to its residues
# (gaps and missing data removed, uppercased, T converted to U)
sub read_stockholm_residues {
  my ($file) = @_;
  my @alisA = ();
  my %seqH  = ();
  my ($line, $name, $seq);
  open(ALI, $file) || return ();
  while($line = <ALI>) {
    if($line =~ /^\/\//) {
      push(@alisA, { %seqH });
      %seqH = ();
    }
    elsif($line =~ /^\#/ || $line !~ /\S/) {
      next;
    }
    elsif($line =~ /^(\S+)\s+(\S+)/) {
      ($name, $seq) = ($1, $2);
      $seq =~ tr/a-z/A-Z/;
      $seq =~ s/[^A-Z]//g;
      $seq =~ tr/T/U/;
      $seqH{$name} .= $seq;
    }
  }
  close(ALI);
  return @alisA;
}
