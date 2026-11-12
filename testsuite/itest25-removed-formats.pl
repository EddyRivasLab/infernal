#! /usr/bin/perl

# Test that the development CM file formats INFERNAL1/b, /c and /d,
# which were never in a release, are rejected cleanly, and that
# cmconvert fails before any work on a CM read from an INFERNAL1/a
# file that would gain something 1/a can't store. (Every other program
# refuses 1/a files: itest26.)
#
# 1. An ASCII file with a 1/b, 1/c or 1/d format tag: every program that
#    reads CMs fails (nonzero exit, no crash) with a message naming the
#    format, saying it was never in a release and to use cmconvert
#    --upgrade.
# 2. The same for binary files with the 1/b, 1/c and 1/d magic numbers,
#    and for a pressed database with them (cmscan, and cmsearch, which
#    reads a CM file's pressed index if it has one), which asks for the
#    file to be pressed again.
# 3. cmconvert --outfmt 1/b and --outfmt 1/d are rejected.
# 4. cmconvert --p7pad of a 1/a file (no --upgrade), --outfmt 1/a
#    --p7pad, and --outfmt 1/e of a 1/a file fail, writing no CM.
#
# Usage:    ./itest25-removed-formats.pl <builddir> <srcdir> <tmpfile prefix>
# Example:  ./itest25-removed-formats.pl ..         ..       tmpfoo
#
# EPN, Fri Nov  6 2026 (w/Claude)
# briefs 26_0821-185/26_0821-186

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
my $lib  = "$srcdir/testsuite/4.1p1.c.cm"; # INFERNAL1/a
my $lib1e = "$srcdir/testsuite/4.c.cm";    # the same CMs, INFERNAL1/e
my $seqs = "$srcdir/testsuite/1k-4.fa";
if (! -e $lib)   { die "FAIL: didn't find $lib\n"; }
if (! -e $lib1e) { die "FAIL: didn't find $lib1e\n"; }
if (! -e $seqs)  { die "FAIL: didn't find $seqs\n"; }
if (first_line($lib)   !~ /^INFERNAL1\/a/) { die "FAIL: $lib is not an INFERNAL1/a file\n"; }
if (first_line($lib1e) !~ /^INFERNAL1\/e/) { die "FAIL: $lib1e is not an INFERNAL1/e file\n"; }
my ($name) = map { /^NAME\s+(\S+)/ ? ($1) : () } read_lines($lib);

my @tmpfiles = ();
sub tmp { my $f = "$tmppfx.$_[0]"; push @tmpfiles, $f; unlink $f if -e $f; return $f; }

my %magic = ("b" => 0xe3edb0b3, "c" => 0xe3edb0b4, "d" => 0xe3edb0b5);
my $err   = tmp("err");

# programs that read a CM file, and how to run each on <cmfile>
sub readers {
    my ($cm) = @_;
    return ("cmalign $cm $seqs",
            "cmcalibrate --forecast $cm",
            "cmconvert $cm",
            "cmemit -N 1 $cm",
            "cmfetch $cm $name",
            "cmpress -F $cm",
            "cmsearch $cm $seqs",
            "cmstat $cm");
}

# 1. ASCII
my $bin1a = tmp("1a.bin");
`$builddir/src/cmconvert -b $lib > $bin1a 2> /dev/null`;
if ($? != 0) { die "FAIL: cmconvert -b failed\n"; }
foreach my $x ("b", "c", "d") {
    my $asc = tmp("1$x.cm");
    write_lines($asc, map { s/^INFERNAL1\/a /INFERNAL1\/$x /r } read_lines($lib));
    foreach my $cmd (readers($asc)) { check_rejected($cmd, $x, "subtest 1"); }

    # 2. binary, and a pressed database
    my $bin = tmp("1$x.bin");
    my $b   = slurp_bin($bin1a);
    substr($b, 0, 4) = pack("L", $magic{$x});
    write_bin($bin, $b);
    foreach my $cmd (readers($bin)) {
        next if $cmd =~ /^cmpress/; # cmpress needs an ASCII file
        check_rejected($cmd, $x, "subtest 2");
    }
    my $db = tmp("db1$x.cm");
    `cp $lib1e $db`;
    foreach my $sfx ("i1m", "i1i", "i1f", "i1p", "i1hf", "i1hp") { tmp("db1$x.cm.$sfx"); }
    `$builddir/src/cmpress -F $db > /dev/null 2>&1`;
    if ($? != 0) { die "FAIL: subtest 2: cmpress of the library failed\n"; }
    $b = slurp_bin("$db.i1m");
    substr($b, 0, 4) = pack("L", $magic{$x});
    write_bin("$db.i1m", $b);
    check_rejected("cmscan $db $seqs",   $x, "subtest 2", 1);
    # a CM file with a dev-format pressed index (as an earlier development cmpress wrote it)
    check_rejected("cmsearch $db $seqs", $x, "subtest 2", 1);
}

# 3.
foreach my $x ("b", "d") {
    `$builddir/src/cmconvert --outfmt 1/$x $lib > /dev/null 2> $err`;
    if ($? == 0) { die "FAIL: subtest 3: cmconvert --outfmt 1/$x didn't fail\n"; }
    if (slurp($err) !~ /No such output format code 1\/$x/) { die "FAIL: subtest 3: cmconvert --outfmt 1/$x: unexpected error:\n" . slurp($err); }
}

# 4.
my $out = tmp("conv.cm");
foreach my $t (["--p7pad --p7pad-N 50",              qr/INFERNAL1\/a.*p7 band pads.*cmconvert --upgrade/s],
               ["--outfmt 1/a --p7pad --p7pad-N 50", qr/--p7pad.*INFERNAL1\/a/s],
               ["--outfmt 1/e",                      qr/--upgrade/]) {
    my ($opts, $re) = @{$t};
    `$builddir/src/cmconvert $opts $lib > $out 2> $err`;
    if ($? == 0)         { die "FAIL: subtest 4: cmconvert $opts of a 1/a file didn't fail\n"; }
    if (($? & 127) != 0) { die "FAIL: subtest 4: cmconvert $opts of a 1/a file crashed\n"; }
    if (slurp($err) !~ $re) { die "FAIL: subtest 4: cmconvert $opts: unexpected error:\n" . slurp($err); }
    if (grep { /^INFERNAL1/ } read_lines($out)) { die "FAIL: subtest 4: cmconvert $opts wrote a CM\n"; }
}

foreach my $f (@tmpfiles) { unlink $f if -e $f; }
print "ok\n";
exit 0;

sub check_rejected {
    my ($cmd, $x, $subtest, $pressed) = @_;
    `$builddir/src/$cmd > $err 2>&1`;
    my $rc = $?;
    my $msg = slurp($err);
    if ($rc == 0)                         { die "FAIL: $subtest: '$cmd' of a 1/$x file didn't fail\n"; }
    if (($rc & 127) != 0)                 { die "FAIL: $subtest: '$cmd' of a 1/$x file crashed\n"; }
    my $ok = $pressed ? ($msg =~ /Pressed file in INFERNAL1\/$x, a never-released development format/ && $msg =~ /cmpress -F/)
                      : ($msg =~ /INFERNAL1\/$x: a development format, never released/ && $msg =~ /cmconvert --upgrade/);
    if (! $ok) {
        die "FAIL: $subtest: '$cmd' of a 1/$x file: unexpected error:\n$msg\n";
    }
}
sub read_lines  { my ($f) = @_; open(my $fh, "<", $f) || die "FAIL: couldn't open $f\n"; my @l = <$fh>; close $fh; return @l; }
sub write_lines { my ($f, @l) = @_; open(my $fh, ">", $f) || die "FAIL: couldn't write $f\n"; print $fh @l; close $fh; }
sub first_line  { my @l = read_lines($_[0]); return (scalar(@l) > 0) ? $l[0] : ""; }
sub slurp       { return join("", read_lines($_[0])); }
sub slurp_bin   { my ($f) = @_; open(my $fh, "<:raw", $f) || die "FAIL: couldn't open $f\n"; local $/; my $b = <$fh>; close $fh; return $b; }
sub write_bin   { my ($f, $b) = @_; open(my $fh, ">:raw", $f) || die "FAIL: couldn't write $f\n"; print $fh $b; close $fh; }
