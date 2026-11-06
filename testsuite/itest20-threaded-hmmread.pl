#! /usr/bin/perl

# Test threaded reads of a pressed INFERNAL1/e CM database, as a
# threaded cmscan does them, on the tutorial's minifam.cm, a released
# INFERNAL1/a file upgraded with cmconvert --upgrade (brief
# 26_0821-186; upgrading it again changes nothing).
#
# A 1/e CM with an HMM-only HMM stores two HMMs. Reading the CM used
# to find the second at the shared HMM file's position after another
# thread could have moved it, so a threaded cmscan of an upgraded
# Rfam library failed with "bad file format for HMM filter".
#
# 1. cmconvert --upgrade on tutorial/minifam.cm (1/e) gives a 1/e file
#    with HMM-only HMMs; cmpress presses it.
# 2. itest_hmmread reads it from several threads (one reading CMs, the
#    others filter HMMs): every read must succeed and give the right
#    HMMs. Before the fix this failed in nearly every run.
# 3. cmscan with threads gives the same hits as without.
#
# Usage:    ./itest20-threaded-hmmread.pl <builddir> <srcdir> <tmpfile prefix>
# Example:  ./itest20-threaded-hmmread.pl ..         ..       tmpfoo
#
# brief 26_0821-168

BEGIN {
    $builddir = shift;
    $srcdir   = shift;
    $tmppfx   = shift;
}

use strict;
use warnings;

our ($builddir, $srcdir, $tmppfx);

for my $prog ("cmconvert", "cmpress", "cmscan", "itest_hmmread") {
    if (! -x "$builddir/src/$prog") { die "FAIL: didn't find $prog binary in $builddir/src\n"; }
}
my $oldcm = "$srcdir/tutorial/minifam.cm";
my $seqs  = "$srcdir/tutorial/metag-example.fa";
if (! -e $oldcm) { die "FAIL: didn't find $oldcm\n"; }
if (! -e $seqs)  { die "FAIL: didn't find $seqs\n"; }

my $cm = "$tmppfx.cm";
unlink glob("$cm.i1*");

# 1. upgrade and press
system("$builddir/src/cmconvert --upgrade --cpu 1 $oldcm > $cm 2> $tmppfx.err");
if ($? != 0) { die "FAIL: cmconvert --upgrade failed\n"; }
open(my $fh, "<", $cm) or die "FAIL: can't open $cm\n";
my $line1 = <$fh>;
my $nh3   = 0;
while (my $line = <$fh>) { if ($line =~ /^HP7\s+h3pri/) { $nh3++; } }
close($fh);
if ($line1 !~ /^INFERNAL1\/e/) { die "FAIL: upgraded file is not in format 1/e\n"; }
if ($nh3 == 0)                  { die "FAIL: upgraded file has no HMM-only HMMs\n"; }
system("$builddir/src/cmpress $cm > /dev/null 2>&1");
if ($? != 0) { die "FAIL: cmpress failed\n"; }

# 2. threaded reads
my $out = `$builddir/src/itest_hmmread -n 4 -N 2000 $cm 2>&1`;
if ($? != 0) { die "FAIL: threaded reads of the pressed database failed:\n$out"; }

# 3. cmscan with and without threads
system("$builddir/src/cmscan --cpu 4 --tblout $tmppfx.t4 $cm $seqs > /dev/null 2>&1");
if ($? != 0) { die "FAIL: cmscan --cpu 4 failed\n"; }
system("$builddir/src/cmscan --cpu 0 --tblout $tmppfx.t0 $cm $seqs > /dev/null 2>&1");
if ($? != 0) { die "FAIL: cmscan --cpu 0 failed\n"; }
my $hits4 = join("", grep { !/^#/ } `cat $tmppfx.t4`);
my $hits0 = join("", grep { !/^#/ } `cat $tmppfx.t0`);
if ($hits4 eq "")     { die "FAIL: cmscan found no hits\n"; }
if ($hits4 ne $hits0) { die "FAIL: cmscan --cpu 4 and --cpu 0 hits differ\n"; }

unlink $cm, glob("$cm.i1*"), "$tmppfx.err", "$tmppfx.t4", "$tmppfx.t0";
print "ok\n";
exit 0;
