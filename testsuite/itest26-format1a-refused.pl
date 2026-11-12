#! /usr/bin/perl

# Test that every program except cmconvert refuses an INFERNAL1/a
# (Infernal 1.1) CM file cleanly, saying to convert it with cmconvert
# --upgrade, as Infernal 1.1 refused v1.0 files outside cmconvert:
#
# 1. An ASCII 1/a file: every program that reads CMs, except cmconvert,
#    fails (nonzero exit, no crash) with that message.
# 2. The same for a binary 1/a file, and for a pressed database in 1/a
#    (cmscan, and cmsearch, which reads a CM file's pressed index if it
#    has one), which asks for the model file to be converted and then
#    pressed again.
# 3. cmconvert reads the ASCII and binary 1/a files.
#
# Usage:    ./itest26-format1a-refused.pl <builddir> <srcdir> <tmpfile prefix>
# Example:  ./itest26-format1a-refused.pl ..         ..       tmpfoo
#
# EPN, Fri Nov  6 2026 (w/Claude)
# brief 26_0821-186

BEGIN {
    $builddir = shift;
    $srcdir   = shift;
    $tmppfx   = shift;
}

use strict;
use warnings;

our ($builddir, $srcdir, $tmppfx);

for my $prog ("cmalign", "cmcalibrate", "cmconvert", "cmemit", "cmfetch", "cmpress", "cmscan", "cmsearch", "cmstat") {
    if (! -x "$builddir/src/$prog") { die "FAIL: didn't find $prog binary in $builddir/src\n"; }
}
my $lib   = "$srcdir/testsuite/4.1p1.c.cm"; # INFERNAL1/a, as Infernal 1.1 wrote it
my $lib1e = "$srcdir/testsuite/4.c.cm";     # the same CMs, INFERNAL1/e
my $seqs  = "$srcdir/testsuite/1k-4.fa";
if (! -e $lib)   { die "FAIL: didn't find $lib\n"; }
if (! -e $lib1e) { die "FAIL: didn't find $lib1e\n"; }
if (! -e $seqs)  { die "FAIL: didn't find $seqs\n"; }
if (first_line($lib)   !~ /^INFERNAL1\/a/) { die "FAIL: $lib is not an INFERNAL1/a file\n"; }
if (first_line($lib1e) !~ /^INFERNAL1\/e/) { die "FAIL: $lib1e is not an INFERNAL1/e file\n"; }
my ($name) = map { /^NAME\s+(\S+)/ ? ($1) : () } read_lines($lib);

my @tmpfiles = ();
sub tmp { my $f = "$tmppfx.$_[0]"; push @tmpfiles, $f; unlink $f if -e $f; return $f; }

my $v1a_magic = 0xe3edb0b2;
my $err       = tmp("err");

# programs other than cmconvert that read a CM file, and how to run each on <cmfile>
sub readers {
    my ($cm) = @_;
    return ("cmalign $cm $seqs",
            "cmcalibrate --forecast $cm",
            "cmcalibrate --nonull3 -L 0.01 $cm",
            "cmemit -N 1 $cm",
            "cmfetch $cm $name",
            "cmpress -F $cm",
            "cmsearch $cm $seqs",
            "cmstat $cm");
}

# 1. ASCII
my $asc = tmp("1a.cm");
`cp $lib $asc`;
foreach my $cmd (readers($asc)) { check_rejected($cmd, "subtest 1"); }
if (slurp($asc) ne slurp($lib)) { die "FAIL: subtest 1: a program changed the 1/a file\n"; }
foreach my $sfx ("i1m", "i1i", "i1f", "i1p", "i1hf", "i1hp", "xxx") {
    if (-e "$asc.$sfx") { die "FAIL: subtest 1: a program wrote $asc.$sfx\n"; }
}

# 2. binary, and a pressed database
my $bin = tmp("1a.bin");
`$builddir/src/cmconvert -b $lib > $bin 2> /dev/null`;
if ($? != 0) { die "FAIL: subtest 2: cmconvert -b of the 1/a file failed\n"; }
if (unpack("L", substr(slurp_bin($bin), 0, 4)) != $v1a_magic) { die "FAIL: subtest 2: cmconvert -b didn't write a binary 1/a file\n"; }
foreach my $cmd (readers($bin)) {
    next if $cmd =~ /^cmpress/; # cmpress needs an ASCII file
    check_rejected($cmd, "subtest 2");
}
# a pressed database whose .i1m file is in 1/a, as Infernal 1.1's cmpress wrote it
my $db = tmp("db1a.cm");
`cp $lib1e $db`;
foreach my $sfx ("i1m", "i1i", "i1f", "i1p", "i1hf", "i1hp") { tmp("db1a.cm.$sfx"); }
`$builddir/src/cmpress -F $db > /dev/null 2>&1`;
if ($? != 0) { die "FAIL: subtest 2: cmpress of the library failed\n"; }
my $b = slurp_bin("$db.i1m");
substr($b, 0, 4) = pack("L", $v1a_magic);
write_bin("$db.i1m", $b);
check_rejected("cmscan $db $seqs",   "subtest 2", 1);
check_rejected("cmsearch $db $seqs", "subtest 2", 1);

# 3.
foreach my $f ($asc, $bin) {
    my $out = tmp("conv.cm");
    `$builddir/src/cmconvert $f > $out 2> $err`;
    if ($? != 0) { die "FAIL: subtest 3: cmconvert of a 1/a file failed:\n" . slurp($err); }
    if (first_line($out) !~ /^INFERNAL1\/a/) { die "FAIL: subtest 3: cmconvert of a 1/a file didn't write INFERNAL1/a\n"; }
    `$builddir/src/cmconvert --upgrade $f > $out 2> $err`;
    if ($? != 0) { die "FAIL: subtest 3: cmconvert --upgrade of a 1/a file failed:\n" . slurp($err); }
    if (first_line($out) !~ /^INFERNAL1\/e/) { die "FAIL: subtest 3: cmconvert --upgrade of a 1/a file didn't write INFERNAL1/e\n"; }
}

foreach my $f (@tmpfiles) { unlink $f if -e $f; }
print "ok\n";
exit 0;

sub check_rejected {
    my ($cmd, $subtest, $pressed) = @_;
    `$builddir/src/$cmd > $err 2>&1`;
    my $rc  = $?;
    my $msg = slurp($err);
    if ($rc == 0)         { die "FAIL: $subtest: '$cmd' of a 1/a file didn't fail\n"; }
    if (($rc & 127) != 0) { die "FAIL: $subtest: '$cmd' of a 1/a file crashed\n"; }
    my $ok = $pressed ? ($msg =~ /INFERNAL1\/a \(Infernal 1\.1\) pressed file: convert the model file with cmconvert --upgrade, then rerun cmpress\./)
                      : ($msg =~ /INFERNAL1\/a \(Infernal 1\.1\) model file: convert it with cmconvert --upgrade first\./);
    if (! $ok) { die "FAIL: $subtest: '$cmd' of a 1/a file: unexpected error:\n$msg\n"; }
}
sub read_lines { my ($f) = @_; open(my $fh, "<", $f) || die "FAIL: couldn't open $f\n"; my @l = <$fh>; close $fh; return @l; }
sub first_line { my @l = read_lines($_[0]); return (scalar(@l) > 0) ? $l[0] : ""; }
sub slurp      { return join("", read_lines($_[0])); }
sub slurp_bin  { my ($f) = @_; open(my $fh, "<:raw", $f) || die "FAIL: couldn't open $f\n"; local $/; my $b = <$fh>; close $fh; return $b; }
sub write_bin  { my ($f, $b) = @_; open(my $fh, ">:raw", $f) || die "FAIL: couldn't write $f\n"; print $fh $b; close $fh; }
