#! /usr/bin/perl

# Test that INFERNAL1/a (Infernal 1.1) CM files are read only by
# cmconvert, which can still write them, for Infernal 1.1:
#
# 1. cmfetch refuses a 1/a library, saying to convert it with cmconvert
#    --upgrade (as every program except cmconvert does: itest26).
# 2. cmconvert (no --outfmt, no --upgrade) of the library writes
#    INFERNAL1/a, equal to the library apart from version stamps,
#    whitespace and the COM line cmconvert adds.
# 3. cmconvert -b, then cmconvert -a, round-trips to the same 1/a file.
# 4. cmconvert --outfmt 1/a of the same CMs in INFERNAL1/e writes
#    INFERNAL1/a, and cmconvert --upgrade of that gives back the 1/e
#    file, apart from version stamps, whitespace and COM lines.
# 5. A 1/a CM with an invalid E-value database size (as some old files
#    have) is refused by cmconvert with a message saying to recalibrate
#    it, and no CM is written.
# (cmconvert of a 1/a CM that would gain something 1/a can't store:
# itest25.)
#
# Usage:    ./itest23-format1a-cmconvert.pl <builddir> <srcdir> <tmpfile prefix>
# Example:  ./itest23-format1a-cmconvert.pl ..         ..       tmpfoo
#
# briefs 26_0821-171/26_0821-185/26_0821-186

BEGIN {
    $builddir = shift;
    $srcdir   = shift;
    $tmppfx   = shift;
}

use strict;
use warnings;

our ($builddir, $srcdir, $tmppfx);

for my $prog ("cmconvert", "cmfetch") {
    if (! -x "$builddir/src/$prog") { die "FAIL: didn't find $prog binary in $builddir/src\n"; }
}
my $lib   = "$srcdir/testsuite/4.1p1.c.cm"; # INFERNAL1/a, as Infernal 1.1 wrote it
my $lib1e = "$srcdir/testsuite/4.c.cm";     # the same CMs, INFERNAL1/e
if (! -e $lib)   { die "FAIL: didn't find $lib\n"; }
if (! -e $lib1e) { die "FAIL: didn't find $lib1e\n"; }
if (first_line($lib)   !~ /^INFERNAL1\/a/) { die "FAIL: $lib is not an INFERNAL1/a file\n"; }
if (first_line($lib1e) !~ /^INFERNAL1\/e/) { die "FAIL: $lib1e is not an INFERNAL1/e file\n"; }

my @tmpfiles = ();
sub tmp { my $f = "$tmppfx.$_[0]"; push @tmpfiles, $f; unlink $f if -e $f; return $f; }
my $err = tmp("err");

# 1.
my ($name) = map { /^NAME\s+(\S+)/ ? ($1) : () } read_lines($lib);
my $one    = tmp("fetch.cm");
`$builddir/src/cmfetch $lib $name > $one 2> $err`;
if ($? == 0)         { die "FAIL: subtest 1: cmfetch of a 1/a library didn't fail\n"; }
if (($? & 127) != 0) { die "FAIL: subtest 1: cmfetch of a 1/a library crashed\n"; }
if (slurp($err) !~ /INFERNAL1\/a \(Infernal 1\.1\) model file: convert it with cmconvert --upgrade first\./) { die "FAIL: subtest 1: unexpected error:\n" . slurp($err); }
if (grep { /^INFERNAL1/ } read_lines($one)) { die "FAIL: subtest 1: cmfetch wrote a CM\n"; }

# 2.
my $conv = tmp("conv.cm");
`$builddir/src/cmconvert $lib > $conv 2> /dev/null`;
if ($? != 0) { die "FAIL: subtest 2: cmconvert failed\n"; }
if (grep { /^INFERNAL1/ && ! /^INFERNAL1\/a/ } read_lines($conv)) { die "FAIL: subtest 2: cmconvert didn't write INFERNAL1/a\n"; }
if (normalized($conv, 1) ne normalized($lib, 1)) { die "FAIL: subtest 2: converted library differs from the 1/a library\n"; }

# 3.
my $bin = tmp("conv.bin");
my $asc = tmp("conv.rt.cm");
`$builddir/src/cmconvert -b $lib > $bin 2> /dev/null`;
if ($? != 0) { die "FAIL: subtest 3: cmconvert -b failed\n"; }
`$builddir/src/cmconvert -a $bin > $asc 2> /dev/null`;
if ($? != 0) { die "FAIL: subtest 3: cmconvert -a failed\n"; }
if (grep { /^INFERNAL1/ && ! /^INFERNAL1\/a/ } read_lines($asc)) { die "FAIL: subtest 3: binary round trip didn't stay INFERNAL1/a\n"; }
if (normalized($asc, 1) ne normalized($lib, 1)) { die "FAIL: subtest 3: binary round trip changed the models\n"; }

# 4.
my $to1a = tmp("to1a.cm");
my $back = tmp("back.cm");
`$builddir/src/cmconvert --outfmt 1/a $lib1e > $to1a 2> $err`;
if ($? != 0) { die "FAIL: subtest 4: cmconvert --outfmt 1/a failed:\n" . slurp($err); }
if (grep { /^INFERNAL1/ && ! /^INFERNAL1\/a/ } read_lines($to1a)) { die "FAIL: subtest 4: cmconvert --outfmt 1/a didn't write INFERNAL1/a\n"; }
`$builddir/src/cmconvert --cpu 1 --upgrade $to1a > $back 2> $err`;
if ($? != 0) { die "FAIL: subtest 4: cmconvert --upgrade of the 1/a file failed:\n" . slurp($err); }
if (normalized($back, 1) ne normalized($lib1e, 1)) { die "FAIL: subtest 4: 1/e -> 1/a -> 1/e changed the models\n"; }

# 5.
my $bad = tmp("bad.cm");
my $badout = tmp("bad.out.cm");
my $nsub = 0;
write_lines($bad, map { ($nsub == 0 && s/^(ECMLC(\s+\S+){3}\s+)\d+/${1}4699193262664056832/) ? do { $nsub++; $_ } : $_ } read_lines($lib));
if ($nsub != 1) { die "FAIL: subtest 5: couldn't make the bad-dbsize file\n"; }
foreach my $opts ("", "--upgrade") {
    `$builddir/src/cmconvert $opts $bad > $badout 2> $err`;
    if ($? == 0)         { die "FAIL: subtest 5: cmconvert $opts of a CM with a bad database size didn't fail\n"; }
    if (($? & 127) != 0) { die "FAIL: subtest 5: cmconvert $opts of a CM with a bad database size crashed\n"; }
    if (slurp($err) !~ /invalid E-value parameters \(database size 4699193262664056832\).*Recalibrate it first/s) { die "FAIL: subtest 5: unexpected error:\n" . slurp($err); }
    if (grep { /^INFERNAL1/ } read_lines($badout)) { die "FAIL: subtest 5: cmconvert $opts wrote a CM\n"; }
}

foreach my $f (@tmpfiles) { unlink $f if -e $f; }
print "ok\n";
exit 0;

sub read_lines { my ($f) = @_; open(my $fh, "<", $f) || die "FAIL: couldn't open $f\n"; my @l = <$fh>; close $fh; return @l; }
sub first_line { my @l = read_lines($_[0]); return (scalar(@l) > 0) ? $l[0] : ""; }
sub slurp      { return join("", read_lines($_[0])); }
sub write_lines { my ($f, @l) = @_; open(my $fh, ">", $f) || die "FAIL: couldn't write $f\n"; print $fh @l; close $fh; }
# file contents without version-stamp lines (and COM lines if <nocom>),
# with each run of whitespace collapsed to one space
sub normalized {
    my ($f, $nocom) = @_;
    my @l = grep { ! /^(INFERNAL1|HMMER3)\// } read_lines($f);
    if ($nocom) { @l = grep { ! /^COM / } @l; }
    return join("", map { s/[ \t]+/ /gr } @l);
}
