#! /usr/bin/perl

# Test that cmbuild and cmconvert --upgrade simulate the glocal
# Forward filter stats (as cmbuild --Esim does) of an entropy-modified
# CM with basepairs whose relative entropy is below the default target,
# and say so:
#
# 1. cmbuild --ere 0.3 on Vault.sto (relative entropy just below the
#    default target) prints one NOTE line, and writes the same CM file
#    as the same build with --Esim (apart from DATE and COM lines).
# 2. A default build prints no such NOTE.
# 3. The CM from 1, written in INFERNAL1/a and converted back with
#    cmconvert --upgrade, gets one NOTE on stderr and the same EFP7GF
#    and EHP7GF lines as cmbuild wrote.
#
# Usage:    ./itest27-lowre-esim.pl <builddir> <srcdir> <tmpfile prefix>
# Example:  ./itest27-lowre-esim.pl ..         ..       tmpfoo
#
# EPN, Tue Nov 24 2026 (w/Claude)
# brief 26_0824-099

BEGIN {
    $builddir = shift;
    $srcdir   = shift;
    $tmppfx   = shift;
}

use strict;
use warnings;

our ($builddir, $srcdir, $tmppfx);

for my $prog ("cmbuild", "cmconvert") {
    if (! -x "$builddir/src/$prog") { die "FAIL: didn't find $prog binary in $builddir/src\n"; }
}
my $aln = "$srcdir/testsuite/Vault.sto";
if (! -e $aln) { die "FAIL: didn't find $aln\n"; }

my @tmpfiles = ();
sub tmp { my $f = "$tmppfx.$_[0]"; push @tmpfiles, $f; unlink $f if -e $f; return $f; }

my $note = qr/^# NOTE: Vault: relative entropy \S+ bits is below the default target \S+ bits; glocal Fwd filter stats simulated/;
my $err  = tmp("err");

# 1.
my ($cm, $out)         = (tmp("lo.cm"),   tmp("lo.out"));
my ($cm_esim, $out_es) = (tmp("es.cm"),   tmp("es.out"));
`$builddir/src/cmbuild --cpu 1 --ere 0.3 -F $cm $aln > $out 2> $err`;
if ($? != 0) { die "FAIL: subtest 1: cmbuild --ere 0.3 failed:\n" . slurp($err); }
my $n = grep { $_ =~ $note } read_lines($out);
if ($n != 1) { die "FAIL: subtest 1: cmbuild --ere 0.3 printed $n NOTE lines, expected 1\n"; }
`$builddir/src/cmbuild --cpu 1 --ere 0.3 --Esim -F $cm_esim $aln > $out_es 2> $err`;
if ($? != 0) { die "FAIL: subtest 1: cmbuild --ere 0.3 --Esim failed:\n" . slurp($err); }
if (strip($cm) ne strip($cm_esim)) { die "FAIL: subtest 1: cmbuild --ere 0.3 and --ere 0.3 --Esim wrote different CMs\n"; }

# 2.
my ($cm_def, $out_def) = (tmp("def.cm"), tmp("def.out"));
`$builddir/src/cmbuild --cpu 1 -F $cm_def $aln > $out_def 2> $err`;
if ($? != 0) { die "FAIL: subtest 2: default cmbuild failed:\n" . slurp($err); }
if (grep { $_ =~ /below the default target/ } read_lines($out_def)) { die "FAIL: subtest 2: default build printed the NOTE\n"; }

# 3.
my ($cm1a, $cmup) = (tmp("1a.cm"), tmp("up.cm"));
`$builddir/src/cmconvert --outfmt 1/a $cm > $cm1a 2> $err`;
if ($? != 0) { die "FAIL: subtest 3: cmconvert --outfmt 1/a failed:\n" . slurp($err); }
`$builddir/src/cmconvert --cpu 1 --upgrade $cm1a > $cmup 2> $err`;
if ($? != 0) { die "FAIL: subtest 3: cmconvert --upgrade failed:\n" . slurp($err); }
$n = grep { $_ =~ $note } read_lines($err);
if ($n != 1) { die "FAIL: subtest 3: cmconvert --upgrade printed $n NOTE lines, expected 1\n"; }
if (hmmstats($cmup) ne hmmstats($cm)) { die "FAIL: subtest 3: cmconvert --upgrade EFP7GF/EHP7GF differ from cmbuild's:\n" . hmmstats($cm) . hmmstats($cmup); }

foreach my $f (@tmpfiles) { unlink $f if -e $f; }
print "ok\n";
exit 0;

sub read_lines { my ($f) = @_; open(my $fh, "<", $f) || die "FAIL: couldn't open $f\n"; my @l = <$fh>; close $fh; return @l; }
sub strip      { return join("", grep { $_ !~ /^(DATE|COM)\s/ } read_lines($_[0])); }
sub hmmstats   { return join("", grep { $_ =~ /^(EFP7GF|EHP7GF)\s/ } read_lines($_[0])); }
