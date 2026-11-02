#! /usr/bin/perl

# Test cmconvert --upgrade, which converts CMs in older formats (with
# a stored filter HMM) to the current format, INFERNAL1/e, deriving and
# calibrating the filter and HMM-only HMMs as cmbuild does; and the
# HP7/HP7A lines that record how the HMM-only HMM is derived.
#
# 1. cmbuild -> cmconvert --outfmt 1/b -> cmconvert --upgrade gives the
#    file cmbuild wrote (apart from DATE and COM lines), for a default
#    build, --p56, --noss (zero basepairs, HP7 ml), and non-default
#    filter HMM calibration options, which --upgrade reads from the
#    cmbuild COM line.
# 2. Without --upgrade, an old-format CM stays in an old format.
# 3. A calibrated INFERNAL1/a file upgrades to INFERNAL1/e, keeping its
#    CM parameters and cmcalibrate E-value parameters; cmsearch runs on
#    it, default and --hmmonly; upgrading the result changes nothing.
# 4. HP7 values written only by development versions (a2, a1) are
#    rejected: formats that were never released aren't supported.
# 5. A v1.0 file converts to INFERNAL1/e without --upgrade, as with it;
#    --cp9ascm writes INFERNAL1/e (no stored filter HMM), with no
#    cmcalibrate E-value parameters. (brief 26_0821-171)
#
# Usage:    ./itest19-cmconvert-upgrade.pl <builddir> <srcdir> <tmpfile prefix>
# Example:  ./itest19-cmconvert-upgrade.pl ..         ..       tmpfoo
#
# brief 26_0821-165

BEGIN {
    $builddir = shift;
    $srcdir   = shift;
    $tmppfx   = shift;
}

use strict;
use warnings;

our ($builddir, $srcdir, $tmppfx);

for my $prog ("cmbuild", "cmconvert", "cmsearch") {
    if (! -x "$builddir/src/$prog") { die "FAIL: didn't find $prog binary in $builddir/src\n"; }
}
my $ali   = "$srcdir/testsuite/tRNA.sto";
my $oldcm = "$srcdir/testsuite/4.c.cm";
my $v1p0cm = "$srcdir/testsuite/4.1p0.cm";
my $seqs  = "$srcdir/testsuite/1k-4.fa";
if (! -e $ali)   { die "FAIL: didn't find $ali\n"; }
if (! -e $oldcm) { die "FAIL: didn't find $oldcm\n"; }
if (! -e $v1p0cm) { die "FAIL: didn't find $v1p0cm\n"; }
if (! -e $seqs)  { die "FAIL: didn't find $seqs\n"; }

my @tmpfiles = ();
sub tmp { my $f = "$tmppfx.$_[0]"; push @tmpfiles, $f; unlink $f if -e $f; return $f; }

# 1.
my @builds = ( [ "default", "",                                      "h3pri" ],
               [ "p56",     "--p56",                                 "h3pri" ],
               [ "noss",    "--noss",                                "ml"    ],
               [ "calopts", "--EmN 20 --EvN 20 --ElfN 20 --Egcmult 1.5 --Eseed 7", "h3pri" ] );
foreach my $b (@builds) {
    my ($tag, $opts, $hp7) = @{$b};
    my $direct = tmp("$tag.cm");
    my $old    = tmp("$tag.1b.cm");
    my $up     = tmp("$tag.up.cm");
    `$builddir/src/cmbuild --cpu 1 -F $opts $direct $ali > /dev/null 2>&1`;
    if ($? != 0) { die "FAIL: subtest 1: cmbuild $opts failed\n"; }
    `$builddir/src/cmconvert --outfmt 1/b $direct > $old 2> /dev/null`;
    if ($? != 0) { die "FAIL: subtest 1: cmconvert --outfmt 1/b ($tag) failed\n"; }
    if (first_line($old) !~ /^INFERNAL1\/b/) { die "FAIL: subtest 1: --outfmt 1/b ($tag) didn't write INFERNAL1/b\n"; }
    `$builddir/src/cmconvert --cpu 1 --upgrade $old > $up 2> /dev/null`;
    if ($? != 0) { die "FAIL: subtest 1: cmconvert --upgrade ($tag) failed\n"; }
    if (first_line($up) !~ /^INFERNAL1\/e/) { die "FAIL: subtest 1: --upgrade ($tag) didn't write INFERNAL1/e\n"; }
    if (stripped($direct) ne stripped($up)) { die "FAIL: subtest 1: upgraded CM ($tag) differs from cmbuild's\n"; }
    my $hp7line = header_value($up, "HP7");
    if (! defined $hp7line || $hp7line ne $hp7) { die "FAIL: subtest 1: ($tag) HP7 is not $hp7\n"; }
    if (defined header_value($up, "HP7A") != ($hp7 eq "h3pri")) { die "FAIL: subtest 1: ($tag) HP7A line should be present iff HP7 is h3pri\n"; }
}

# 2.
my $plain = tmp("4.plain.cm");
`$builddir/src/cmconvert $oldcm > $plain 2> /dev/null`;
if ($? != 0) { die "FAIL: subtest 2: cmconvert failed\n"; }
if (first_line($plain) =~ /^INFERNAL1\/e/) { die "FAIL: subtest 2: cmconvert without --upgrade wrote INFERNAL1/e\n"; }

# 3.
my $up4  = tmp("4.up.cm");
my $up4b = tmp("4.up2.cm");
`$builddir/src/cmconvert --cpu 1 --upgrade $oldcm > $up4 2> /dev/null`;
if ($? != 0) { die "FAIL: subtest 3: cmconvert --upgrade failed\n"; }
my @fmt = grep { /^INFERNAL1/ } read_lines($up4);
if (scalar(@fmt) != 4 || grep { ! /^INFERNAL1\/e/ } @fmt) { die "FAIL: subtest 3: expected 4 INFERNAL1/e CMs\n"; }
if (join("", grep { /^ECM/ } read_lines($oldcm)) ne join("", grep { /^ECM/ } read_lines($up4))) { die "FAIL: subtest 3: ECM lines changed\n"; }
if (join("", state_lines($oldcm)) ne join("", state_lines($up4))) { die "FAIL: subtest 3: CM parameters changed\n"; }
foreach my $opts ("", "--hmmonly") {
    `$builddir/src/cmsearch --cpu 0 $opts $up4 $seqs > /dev/null 2>&1`;
    if ($? != 0) { die "FAIL: subtest 3: cmsearch $opts on the upgraded CM failed\n"; }
}
`$builddir/src/cmconvert --upgrade $up4 > $up4b 2> /dev/null`;
if ($? != 0) { die "FAIL: subtest 3: cmconvert --upgrade of a current-format CM failed\n"; }
if (stripped($up4) ne stripped($up4b)) { die "FAIL: subtest 3: upgrading a current-format CM changed it\n"; }

# 4.
my $direct = "$tmppfx.default.cm";
my $a2     = tmp("a2.cm");
my $a1     = tmp("a1.cm");
write_lines($a2, map { s/^HP7(\s+)h3pri/HP7$1a2/r } read_lines($direct));
`$builddir/src/cmconvert $a2 > /dev/null 2>&1`;
if ($? == 0) { die "FAIL: subtest 4: HP7 a2 file not rejected\n"; }
write_lines($a1, map { s/^HP7(\s+)h3pri/HP7$1a1/r } grep { ! /^HP7A / } read_lines($direct));
`$builddir/src/cmconvert $a1 > /dev/null 2>&1`;
if ($? == 0) { die "FAIL: subtest 4: HP7 a1 file not rejected\n"; }

# 5.
my $v1     = tmp("v1p0.cm");
my $v1up   = tmp("v1p0.up.cm");
my $cp9    = tmp("cp9.cm");
`$builddir/src/cmconvert $v1p0cm > $v1 2> /dev/null`;
if ($? != 0) { die "FAIL: subtest 5: cmconvert of a v1.0 file failed\n"; }
`$builddir/src/cmconvert --upgrade $v1p0cm > $v1up 2> /dev/null`;
if ($? != 0) { die "FAIL: subtest 5: cmconvert --upgrade of a v1.0 file failed\n"; }
if (grep { /^INFERNAL1/ && ! /^INFERNAL1\/e/ } read_lines($v1)) { die "FAIL: subtest 5: v1.0 file not converted to INFERNAL1/e\n"; }
if (stripped($v1) ne stripped($v1up)) { die "FAIL: subtest 5: v1.0 conversion differs from --upgrade\n"; }
`$builddir/src/cmconvert --cp9ascm $oldcm > $cp9 2> /dev/null`;
if ($? != 0) { die "FAIL: subtest 5: cmconvert --cp9ascm failed\n"; }
if (grep { /^INFERNAL1/ && ! /^INFERNAL1\/e/ } read_lines($cp9)) { die "FAIL: subtest 5: --cp9ascm didn't write INFERNAL1/e\n"; }
if (grep { /^ECM/ } read_lines($cp9)) { die "FAIL: subtest 5: --cp9ascm CM has cmcalibrate E-value parameters\n"; }
`$builddir/src/cmsearch --cpu 0 --hmmonly $cp9 $seqs > /dev/null 2>&1`;
if ($? != 0) { die "FAIL: subtest 5: cmsearch --hmmonly on the --cp9ascm CMs failed\n"; }

foreach my $f (@tmpfiles) { unlink $f if -e $f; }
print "ok\n";
exit 0;

sub read_lines  { my ($f) = @_; open(my $fh, "<", $f) || die "FAIL: couldn't open $f\n"; my @l = <$fh>; close $fh; return @l; }
sub write_lines { my ($f, @l) = @_; open(my $fh, ">", $f) || die "FAIL: couldn't open $f\n"; print $fh @l; close $fh; }
sub first_line  { my @l = read_lines($_[0]); return (scalar(@l) > 0) ? $l[0] : ""; }
sub stripped    { return join("", grep { ! /^(DATE|COM) / } read_lines($_[0])); }
sub header_value {
    my ($f, $tag) = @_;
    foreach my $l (read_lines($f)) { if ($l =~ /^\Q$tag\E\s+(.*?)\s*$/) { return $1; } if ($l =~ /^CM\s*$/) { last; } }
    return undef;
}
# CM state lines (not node lines, whose trailing columns depend on the format)
sub state_lines {
    my ($f) = @_;
    my $in = 0; my @s = ();
    foreach my $l (read_lines($f)) {
        if ($l =~ /^CM\s*$/) { $in = 1; next; }
        if ($l =~ /^\/\//)   { $in = 0; next; }
        if ($in && $l !~ /\[ /) { push @s, $l; }
    }
    return @s;
}
