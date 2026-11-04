#! /usr/bin/perl

# itest22: cmscan and cmsearch must search each model with the same HMM.
#
# A model that is searched HMM-only (every model with --hmmonly; by
# default, a model with 0 basepairs) is searched with its HMM-only
# search HMM (cm->hp7, if it has one) and that HMM's E-value
# parameters. cmsearch has always done that; cmscan did it only with
# --hmmonly, and searched a 0 basepair model in its default mode with
# its filter HMM instead, which gave different scores, E-values and
# hits (brief 26_0821-177).
#
# We press a mixed database of three models built from tRNA.sto:
#   trnaP: 0 basepairs (--noss), default build (no separate hp7)
#   trnaN: 0 basepairs (--noss --noh3pri), which has an hp7
#   trnaS: with its structure
# and check that for each model, in default, --hmmonly, -g and
# --hmmonly -g mode, cmscan (serial and threaded) reports the same
# hits as cmsearch on that model alone: the same model and sequence
# coordinates, strand, truncation, pass, GC, bias, score, E-value and
# inclusion. -Z is set explicitly, so E-values are comparable. Every
# model must have at least one hit in default mode.
#
# EPN, Tue Nov  3 2026 (w/Claude)
#

# Usage:    ./itest22-scan-hmmonly-hmm.pl <builddir> <srcdir> <tmpfile prefix>
# Example:  ./itest22-scan-hmmonly-hmm.pl ..         ..       tmpfoo

$usage = "perl itest22-scan-hmmonly-hmm.pl <builddir> <srcdir> <tmpfile prefix>\n";
if ($#ARGV != 2) { die "Wrong argument number.\n$usage"; }

$builddir = shift;
$srcdir   = shift;
$tmp      = shift;
$ok       = 1;

foreach $prog ("cmbuild", "cmpress", "cmsearch", "cmscan") { 
  if (! -x "$builddir/src/$prog") { die "FAIL: didn't find $prog binary in $builddir/src\n"; }
}
$cmbuild  = "$builddir/src/cmbuild";
$cmpress  = "$builddir/src/cmpress";
$cmsearch = "$builddir/src/cmsearch";
$cmscan   = "$builddir/src/cmscan";
$alifile  = "$srcdir/testsuite/tRNA.sto";
$seqfile  = "$srcdir/testsuite/emitted-tRNA.fa";
if (! -e $alifile) { die "FAIL: didn't find $alifile\n"; }
if (! -e $seqfile) { die "FAIL: didn't find $seqfile\n"; }

@models = ("trnaP", "trnaN", "trnaS");
%bopts  = ("trnaP" => "--noss", "trnaN" => "--noss --noh3pri", "trnaS" => "");
@modes  = ("", "--hmmonly", "-g", "--hmmonly -g");

# build the models and the pressed database
if ($ok) {
  unlink("$tmp.db.cm");
  foreach $m (@models) {
    `$cmbuild -F -n $m $bopts{$m} $tmp.$m.raw.cm $alifile`;
    if ($? != 0) { $ok = 0; last; }
    # drop the accession (the same for all three models), so cmpress can index them
    `grep -v '^ACC ' $tmp.$m.raw.cm > $tmp.$m.cm`;
    `cat $tmp.$m.cm >> $tmp.db.cm`;
  }
}
if ($ok) {
  `$cmpress -F $tmp.db.cm`;
  if ($? != 0) { $ok = 0; }
}

# compare per model and mode
foreach $mode (@modes) {
  last if (! $ok);
  %scanH = ();
  foreach $cpu ("0", "2") {
    `$cmscan --cpu $cpu -Z 1 $mode --tblout $tmp.scan.tbl $tmp.db.cm $seqfile`;
    if ($? != 0) { $ok = 0; last; }
    %{$scanH{$cpu}} = read_tblout("$tmp.scan.tbl", 1);
  }
  last if (! $ok);
  foreach $m (@models) {
    `$cmsearch --cpu 0 -Z 1 $mode --tblout $tmp.search.tbl $tmp.$m.cm $seqfile`;
    if ($? != 0) { $ok = 0; last; }
    %searchH = read_tblout("$tmp.search.tbl", 0);
    $search  = (exists $searchH{$m}) ? $searchH{$m} : "";
    if ($mode eq "" && $search eq "") {
      printf("no hits to model $m in default mode\n");
      $ok = 0;
    }
    foreach $cpu ("0", "2") {
      $scan = (exists $scanH{$cpu}{$m}) ? $scanH{$cpu}{$m} : "";
      if ($scan ne $search) {
        printf("model $m, mode '$mode', cmscan --cpu $cpu differs from cmsearch:\ncmsearch:\n%scmscan:\n%s", $search, $scan);
        $ok = 0;
      }
    }
  }
}

foreach $sfx ("db.cm", "db.cm.i1m", "db.cm.i1i", "db.cm.i1f", "db.cm.i1p", "db.cm.i1hf", "db.cm.i1hp", "scan.tbl", "search.tbl", "trnaP.cm", "trnaN.cm", "trnaS.cm", "trnaP.raw.cm", "trnaN.raw.cm", "trnaS.raw.cm") {
  unlink "$tmp.$sfx";
}

if ($ok) { print "ok\n"; exit 0; }
print "FAILED\n";
exit 1;

# read_tblout(<file>, <is_scan>): return a hash, key: model name, value:
# that model's hit lines, sorted, with the model and sequence name and
# accession columns in cmsearch order and the description dropped
# (cmscan prints the model's, cmsearch the sequence's).
sub read_tblout {
  my ($file, $is_scan) = @_;
  my %linesH = ();
  my %retH   = ();
  my @el;
  open(IN, $file) || die "ERROR unable to open $file";
  while(<IN>) {
    next if m/^\#/;
    chomp;
    @el = split(/\s+/);
    if ($is_scan) { @el[0,1,2,3] = @el[2,3,0,1]; }
    # now: el[0..1] sequence name/acc, el[2..3] model name/acc
    push(@{$linesH{$el[2]}}, join(" ", @el[0..16]) . "\n");
  }
  close(IN);
  foreach my $m (keys %linesH) {
    $retH{$m} = join("", sort @{$linesH{$m}});
  }
  return %retH;
}
