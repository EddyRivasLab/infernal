#! /usr/bin/perl

# bug i84 (GitHub #84) - --noF4 did not actually turn off the stage 4
# glocal HMM Forward filter
#
# In pli_p7_env_def() (src/cm_pipeline.c), the glocal Forward P-value
# test (stage 4, F4) was applied unconditionally: --noF4 set
# pli->do_gfwd = FALSE, but that flag was only read when printing the
# pipeline statistics, never to skip the actual P > pli->F4 test. A
# search run with --noF4 therefore kept filtering on F4 at its default
# threshold instead of turning the filter off, and could miss hits
# that the same search with the filter genuinely opened (--F4 1.0)
# would find.
#
# bug-i84.cm is the tutorial's tRNA5.c.cm; bug-i84.fa is a 500-residue
# slice of the tutorial's mrum-genome.fa spanning a tRNA hit that
# --noF4's bug filtered out. We check that 'cmsearch --noF1 --noF3
# --noF3b --noF4 --noF4b --F5 1.0 --noF6 --notrunc' (F4 off) finds the
# same hits as the same command with '--F4 1.0 --noF4b' in place of
# '--noF4 --noF4b' (F4 on but open at threshold 1.0 so it filters
# nothing) -- the two must be equivalent once --noF4 genuinely means
# "off".
#
# EPN, Fri Nov 20 2026 (w/Claude)
#

$usage = "perl bug-i84.pl <cmsearch> <testsuite/bug-i84.cm> <testsuite/bug-i84.fa>\n";
if ($#ARGV != 2) { die "Wrong argument number.\n$usage"; }

$cmsearch = shift;
$cmfile   = shift;
$seqfile  = shift;
$ok       = 1;

$common = "--noF1 --noF3 --noF3b --F5 1.0 --noF6 --notrunc";

if ($ok) {
  $output = `$cmsearch $common --noF4 --noF4b --tblout i84.a.tbl $cmfile $seqfile 2>&1`;
  if ($? != 0 || $output !~ /\n\# CPU/) { $ok = 0; }
}

if ($ok) {
  $output = `$cmsearch $common --F4 1.0 --noF4b --tblout i84.b.tbl $cmfile $seqfile 2>&1`;
  if ($? != 0 || $output !~ /\n\# CPU/) { $ok = 0; }
}

if ($ok) {
  @a = read_tblout_hits("i84.a.tbl");
  @b = read_tblout_hits("i84.b.tbl");
  if (scalar(@a) == 0)             { $ok = 0; } # sanity: the region has >= 1 hit either way
  if (scalar(@a) != scalar(@b))    { $ok = 0; } # --noF4 must match --F4 1.0 (both effectively "off")
  else {
    for ($i = 0; $i < scalar(@a); $i++) {
      if ($a[$i] ne $b[$i]) { $ok = 0; }
    }
  }
}

foreach $tmpfile ("i84.a.tbl", "i84.b.tbl") {
  unlink $tmpfile if -e $tmpfile;
}
if ($ok) { print "ok\n";     exit 0; }
else     { print "FAILED\n"; exit 1; }

# read_tblout_hits(<file>): return a sorted list of "seqfrom-seqto" hit
# coordinate strings from a --tblout file (ignoring comment lines)
sub read_tblout_hits {
  my ($file) = @_;
  my @hits = ();
  my ($line);
  open(TBL, $file) || return ();
  while ($line = <TBL>) {
    next if ($line =~ /^\#/);
    next if ($line !~ /\S/);
    my @el = split(/\s+/, $line);
    push(@hits, "$el[7]-$el[8]"); # seq from, seq to columns
  }
  close(TBL);
  return sort(@hits);
}
