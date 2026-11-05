#! /usr/bin/perl

# Test that a CM read from a released INFERNAL1/a file is written back
# as INFERNAL1/a by the programs that write the default format, so that
# Infernal 1.1 can still read it:
#
# 1. cmfetch of each model in a 1/a library writes INFERNAL1/a, and
#    the fetched models together equal the library, apart from the
#    INFERNAL1/HMMER3 version-stamp lines and whitespace.
# 2. cmconvert (no --outfmt, no --upgrade) of the library writes
#    INFERNAL1/a, equal to the library apart from version stamps,
#    whitespace and the COM line cmconvert adds.
# 3. cmconvert -b, then cmconvert -a, round-trips to the same 1/a file.
# 4. A 1/b file stays 1/b under cmfetch and cmconvert.
#
# Usage:    ./itest23-format1a-stays-1a.pl <builddir> <srcdir> <tmpfile prefix>
# Example:  ./itest23-format1a-stays-1a.pl ..         ..       tmpfoo
#
# brief 26_0821-171

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
my $lib = "$srcdir/testsuite/4.c.cm";
if (! -e $lib) { die "FAIL: didn't find $lib\n"; }
if (first_line($lib) !~ /^INFERNAL1\/a/) { die "FAIL: $lib is not an INFERNAL1/a file\n"; }

my @tmpfiles = ();
sub tmp { my $f = "$tmppfx.$_[0]"; push @tmpfiles, $f; unlink $f if -e $f; return $f; }

# 1.
my %seen    = ();
my @names   = grep { ! $seen{$_}++ } map { /^NAME\s+(\S+)/ ? ($1) : () } read_lines($lib);  # CM and filter HMM share a NAME
my $fetched = "";
foreach my $name (@names) {
    my $one = tmp("fetch.cm");
    `$builddir/src/cmfetch $lib $name > $one 2> /dev/null`;
    if ($? != 0) { die "FAIL: subtest 1: cmfetch $name failed\n"; }
    if (first_line($one) !~ /^INFERNAL1\/a/) { die "FAIL: subtest 1: cmfetch $name didn't write INFERNAL1/a\n"; }
    $fetched .= normalized($one);
}
if ($fetched ne normalized($lib)) { die "FAIL: subtest 1: fetched models differ from the 1/a library\n"; }

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
my $old  = tmp("1b.cm");
my $old2 = tmp("1b.conv.cm");
my $old3 = tmp("1b.fetch.cm");
`$builddir/src/cmconvert --outfmt 1/b $lib > $old 2> /dev/null`;
if ($? != 0) { die "FAIL: subtest 4: cmconvert --outfmt 1/b failed\n"; }
`$builddir/src/cmconvert $old > $old2 2> /dev/null`;
if ($? != 0) { die "FAIL: subtest 4: cmconvert of a 1/b file failed\n"; }
if (grep { /^INFERNAL1/ && ! /^INFERNAL1\/b/ } read_lines($old2)) { die "FAIL: subtest 4: cmconvert didn't keep INFERNAL1/b\n"; }
`$builddir/src/cmfetch $old $names[0] > $old3 2> /dev/null`;
if ($? != 0) { die "FAIL: subtest 4: cmfetch of a 1/b file failed\n"; }
if (first_line($old3) !~ /^INFERNAL1\/b/) { die "FAIL: subtest 4: cmfetch didn't keep INFERNAL1/b\n"; }

foreach my $f (@tmpfiles) { unlink $f if -e $f; }
print "ok\n";
exit 0;

sub read_lines { my ($f) = @_; open(my $fh, "<", $f) || die "FAIL: couldn't open $f\n"; my @l = <$fh>; close $fh; return @l; }
sub first_line { my @l = read_lines($_[0]); return (scalar(@l) > 0) ? $l[0] : ""; }
# file contents without version-stamp lines (and COM lines if <nocom>),
# with each run of whitespace collapsed to one space
sub normalized {
    my ($f, $nocom) = @_;
    my @l = grep { ! /^(INFERNAL1|HMMER3)\// } read_lines($f);
    if ($nocom) { @l = grep { ! /^COM / } @l; }
    return join("", map { s/[ \t]+/ /gr } @l);
}
