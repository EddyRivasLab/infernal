#! /usr/bin/perl

# Test the INFERNAL1/d (calibrated cmbuild) and INFERNAL1/b (cmconvert
# of legacy v1a) CM file formats and --p7pad-* options of cmbuild
# and cmconvert. Regression-covers the c09285c8 pass-2 rewind dispatch
# fix (cmcalibrate rereads v1b output) via a separate sqc entry in
# dev_testsuite.sqc at level 2; this script covers the ASCII-format
# invariants, determinism, serial-vs-threaded equivalence, and the
# cmconvert v1a -> v1b pad injection path.
#
# Usage:    ./itest12-p7nodepad.pl <builddir> <srcdir> <tmpfile prefix>
# Example:  ./itest12-p7nodepad.pl ..         ..       tmpfoo
#
# All cmbuild/cmconvert invocations use --p7pad-N 50 (not the default
# 1000) so the Monte Carlo pad simulation stays fast.

BEGIN {
    $builddir = shift;
    $srcdir   = shift;
    $tmppfx   = shift;
}

use strict;
use warnings;

our ($builddir, $srcdir, $tmppfx);

# Verify that we have all the executables and input files we need.
if (! -x "$builddir/src/cmbuild")            { die "FAIL: didn't find cmbuild binary in $builddir/src\n";   }
if (! -x "$builddir/src/cmconvert")          { die "FAIL: didn't find cmconvert binary in $builddir/src\n"; }
if (! -e "$srcdir/testsuite/tRNA.sto")       { die "FAIL: didn't find tRNA.sto in $srcdir/testsuite\n";     }
if (! -e "$srcdir/testsuite/tRNA.1p0.cm")    { die "FAIL: didn't find tRNA.1p0.cm in $srcdir/testsuite\n";  }

my $ali    = "$srcdir/testsuite/tRNA.sto";
my $v1a_in = "$srcdir/testsuite/tRNA.1p0.cm";

# Clean up any leftover files from a previous run
for my $suffix ("default.cm", "nopad.cm",
                "seed17a.cm", "seed17b.cm",
                "serial.cm", "threaded.cm",
                "bin.cm", "back.cm",
                "1p0-input.cm", "converted.cm", "converted.nopad.cm",
                "seed17a.filt", "seed17b.filt",
                "cmbuild.log") {
    if (-e "$tmppfx.$suffix") { unlink "$tmppfx.$suffix"; }
}

######################################################################
# Subtest 1: Default cmbuild (with pad) produces INFERNAL1/d with pads.
# (Default cmbuild now stores a null3-OFF E-value set via the fastcal
# store-both path, so the CM carries CMH_EXPTAIL_NONULL3_STATS and is
# written as INFERNAL1/d, not 1/b -- see cm_file.c write-format gate.)
######################################################################
`$builddir/src/cmbuild --wnone -F --p7pad-N 50 --p7pad-seed 42 $tmppfx.default.cm $ali > $tmppfx.cmbuild.log 2>&1`;
if ($? != 0) { die "FAIL: cmbuild (default) failed\n"; }

my $first = first_line("$tmppfx.default.cm");
if ($first !~ /^INFERNAL1\/d\b/) { die "FAIL: subtest 1: default cmbuild first line is '$first', expected to start with 'INFERNAL1/d'\n"; }
if (! file_has_header("$tmppfx.default.cm", "P7NODEPAD", "yes")) { die "FAIL: subtest 1: default cmbuild output missing 'P7NODEPAD yes'\n"; }

######################################################################
# Subtest 2: --no-p7pad suppresses pad computation.
# First line is INFERNAL1/d (format tag is independent of whether pads
# are populated; default cmbuild still stores the null3-OFF E-value set,
# so the CM is 1/d), header says 'P7NODEPAD no', and the first MATL node
# line trails with '- -' (both dashes) rather than '<int> -'.
######################################################################
`$builddir/src/cmbuild --wnone -F --no-p7pad $tmppfx.nopad.cm $ali > $tmppfx.cmbuild.log 2>&1`;
if ($? != 0) { die "FAIL: cmbuild --no-p7pad failed\n"; }

$first = first_line("$tmppfx.nopad.cm");
if ($first !~ /^INFERNAL1\/d\b/) { die "FAIL: subtest 2: --no-p7pad first line is '$first', expected to start with 'INFERNAL1/d'\n"; }
if (! file_has_header("$tmppfx.nopad.cm", "P7NODEPAD", "no")) { die "FAIL: subtest 2: --no-p7pad output missing 'P7NODEPAD no'\n"; }

# Inspect the first MATL node line: last two whitespace fields must be '- -'.
my $matl_nopad = first_matl_line("$tmppfx.nopad.cm");
if (! defined $matl_nopad) { die "FAIL: subtest 2: no MATL node line found in $tmppfx.nopad.cm\n"; }
my @f_nopad = split /\s+/, $matl_nopad;
# leading empty element from leading whitespace?
shift @f_nopad if $f_nopad[0] eq "";
# INFERNAL1/d node lines carry 14 whitespace fields: the 12 base 1/b
# columns plus two trailing consensus-pseudoknot annotation columns
# (idx 12,13). The left/right p7 pad fields remain at idx 10,11.
if (scalar(@f_nopad) != 14) { die "FAIL: subtest 2: MATL line should have 14 whitespace fields, got " . scalar(@f_nopad) . ": '$matl_nopad'\n"; }
if ($f_nopad[10] ne "-" || $f_nopad[11] ne "-") {
    die "FAIL: subtest 2: --no-p7pad MATL last two fields should be '- -', got '$f_nopad[10]' '$f_nopad[11]'\n";
}

######################################################################
# Subtest 3: With pads, MATL node lines have 14 fields (12 base 1/b
# columns + 2 trailing 1/d pknot columns); for MATL (left-only match)
# the pad fields at idx 10,11 are <int> followed by '-'.
######################################################################
my $matl_pad = first_matl_line("$tmppfx.default.cm");
if (! defined $matl_pad) { die "FAIL: subtest 3: no MATL node line found in $tmppfx.default.cm\n"; }
my @f_pad = split /\s+/, $matl_pad;
shift @f_pad if $f_pad[0] eq "";
if (scalar(@f_pad) != 14) { die "FAIL: subtest 3: MATL line should have 14 whitespace fields, got " . scalar(@f_pad) . ": '$matl_pad'\n"; }
if ($f_pad[10] !~ /^\d+$/) { die "FAIL: subtest 3: MATL 11th field (left pad) should be integer, got '$f_pad[10]'\n"; }
if ($f_pad[11] ne "-")     { die "FAIL: subtest 3: MATL 12th field (right pad) should be '-' for MATL, got '$f_pad[11]'\n"; }

######################################################################
# Subtest 4: Determinism. Two builds with the same --p7pad-N/--p7pad-seed
# must produce identical ASCII CMs (modulo COM lines that embed cmdline
# and timestamps).
#
# We fix --Eseed (HMM filter calibration RNG), --cpu 1 (so the glocal
# Forward tau fit is not subject to thread scheduling), --p7pad-seed,
# and --p7pad-cpu 0 (serial pad Monte Carlo). The threaded path is
# exercised separately as subtest 5.
######################################################################
`$builddir/src/cmbuild --wnone -F --cpu 1 --Eseed 42 --p7pad-N 50 --p7pad-seed 17 --p7pad-cpu 0 $tmppfx.seed17a.cm $ali > $tmppfx.cmbuild.log 2>&1`;
if ($? != 0) { die "FAIL: cmbuild (seed17 a) failed\n"; }
`$builddir/src/cmbuild --wnone -F --cpu 1 --Eseed 42 --p7pad-N 50 --p7pad-seed 17 --p7pad-cpu 0 $tmppfx.seed17b.cm $ali > $tmppfx.cmbuild.log 2>&1`;
if ($? != 0) { die "FAIL: cmbuild (seed17 b) failed\n"; }

strip_com_lines("$tmppfx.seed17a.cm", "$tmppfx.seed17a.filt");
strip_com_lines("$tmppfx.seed17b.cm", "$tmppfx.seed17b.filt");

my $diff = `diff $tmppfx.seed17a.filt $tmppfx.seed17b.filt`;
if ($diff ne "") { die "FAIL: subtest 4: two builds with same seed differ after stripping COM lines:\n$diff\n"; }

######################################################################
# Subtest 5: Serial and threaded pad computation both complete and
# write valid pad fields. Serial-vs-threaded byte-identity is NOT
# guaranteed -- the threaded path seeds one RNG per worker from the
# master RNG and workers pull samples from a queue, so the
# decomposition of the Monte Carlo sample stream differs. Similarly,
# two threaded runs may differ due to queue scheduling order even at
# the same seed and ncpu. So this subtest verifies (a) both paths
# succeed, (b) both produce per-node pad fields that parse as
# integers where required, and (c) the total number of pad rows
# matches between serial and threaded outputs (same model topology).
######################################################################
`$builddir/src/cmbuild --wnone -F --cpu 1 --Eseed 42 --p7pad-N 50 --p7pad-seed 17 --p7pad-cpu 0 $tmppfx.serial.cm   $ali > $tmppfx.cmbuild.log 2>&1`;
if ($? != 0) { die "FAIL: cmbuild (serial) failed\n"; }
`$builddir/src/cmbuild --wnone -F --cpu 1 --Eseed 42 --p7pad-N 50 --p7pad-seed 17 --p7pad-cpu 4 $tmppfx.threaded.cm $ali > $tmppfx.cmbuild.log 2>&1`;
if ($? != 0) { die "FAIL: cmbuild (threaded) failed\n"; }

my $pads_serial   = extract_pad_columns("$tmppfx.serial.cm");
my $pads_threaded = extract_pad_columns("$tmppfx.threaded.cm");
if ($pads_serial eq "")   { die "FAIL: subtest 5: no pad columns extracted from serial build\n"; }
if ($pads_threaded eq "") { die "FAIL: subtest 5: no pad columns extracted from threaded build\n"; }
my $n_serial_rows   = ($pads_serial   =~ tr/\n//);
my $n_threaded_rows = ($pads_threaded =~ tr/\n//);
if ($n_serial_rows != $n_threaded_rows) {
    die "FAIL: subtest 5: serial has $n_serial_rows MAT* rows, threaded has $n_threaded_rows (topology mismatch)\n";
}
# Also confirm both serial and threaded CMs advertise 'P7NODEPAD yes'.
if (! file_has_header("$tmppfx.serial.cm",   "P7NODEPAD", "yes")) { die "FAIL: subtest 5: serial CM missing 'P7NODEPAD yes'\n"; }
if (! file_has_header("$tmppfx.threaded.cm", "P7NODEPAD", "yes")) { die "FAIL: subtest 5: threaded CM missing 'P7NODEPAD yes'\n"; }

######################################################################
# Subtest 6: Binary round-trip preserves pads.
# cmconvert -b writes a binary v1b, then cmconvert -a converts back to
# ASCII. Pad columns on MAT* node lines must match the original.
# We pass --no-p7pad on both conversions to suppress pad re-computation
# (the input already has pads from cmbuild; cmconvert's compute path
# requires a configure_model() pass which is only run for v1.0 inputs,
# so without --no-p7pad this would error out complaining that cp9map
# is missing). With --no-p7pad, the already-embedded pads are read and
# re-written as-is.
######################################################################
`$builddir/src/cmconvert --no-p7pad -b $tmppfx.default.cm > $tmppfx.bin.cm 2>>$tmppfx.cmbuild.log`;
if ($? != 0) { die "FAIL: cmconvert -b (binary round-trip write) failed\n"; }
`$builddir/src/cmconvert --no-p7pad -a $tmppfx.bin.cm > $tmppfx.back.cm 2>>$tmppfx.cmbuild.log`;
if ($? != 0) { die "FAIL: cmconvert -a (binary round-trip read) failed\n"; }

my $pads_orig = extract_pad_columns("$tmppfx.default.cm");
my $pads_back = extract_pad_columns("$tmppfx.back.cm");
if ($pads_orig ne $pads_back) {
    die "FAIL: subtest 6: binary round-trip altered pad columns\n--- original ---\n$pads_orig\n--- round-tripped ---\n$pads_back\n";
}

######################################################################
# Subtest 7: cmconvert v1a -> v1b adds pads.
# Copy the checked-in legacy 1.0 CM locally first because cmconvert
# writes to the current directory and we want isolation.
######################################################################
`cp $v1a_in $tmppfx.1p0-input.cm`;
if ($? != 0) { die "FAIL: failed to copy legacy v1a input CM\n"; }

`$builddir/src/cmconvert --p7pad-N 50 --p7pad-seed 42 $tmppfx.1p0-input.cm > $tmppfx.converted.cm 2>>$tmppfx.cmbuild.log`;
if ($? != 0) { die "FAIL: cmconvert v1a->v1b (with pads) failed\n"; }

$first = first_line("$tmppfx.converted.cm");
if ($first !~ /^INFERNAL1\/b\b/) { die "FAIL: subtest 7: cmconvert output first line is '$first', expected to start with 'INFERNAL1/b'\n"; }
if (! file_has_header("$tmppfx.converted.cm", "P7NODEPAD", "yes")) { die "FAIL: subtest 7: cmconvert output missing 'P7NODEPAD yes'\n"; }

######################################################################
# Subtest 8: cmconvert --no-p7pad produces v1b with pads disabled.
######################################################################
`$builddir/src/cmconvert --no-p7pad $tmppfx.1p0-input.cm > $tmppfx.converted.nopad.cm 2>>$tmppfx.cmbuild.log`;
if ($? != 0) { die "FAIL: cmconvert --no-p7pad failed\n"; }

$first = first_line("$tmppfx.converted.nopad.cm");
if ($first !~ /^INFERNAL1\/b\b/) { die "FAIL: subtest 8: cmconvert --no-p7pad first line is '$first', expected to start with 'INFERNAL1/b'\n"; }
if (! file_has_header("$tmppfx.converted.nopad.cm", "P7NODEPAD", "no")) { die "FAIL: subtest 8: cmconvert --no-p7pad output missing 'P7NODEPAD no'\n"; }

######################################################################

print "ok\n";

# cleanup
for my $suffix ("default.cm", "nopad.cm",
                "seed17a.cm", "seed17b.cm",
                "serial.cm", "threaded.cm",
                "bin.cm", "back.cm",
                "1p0-input.cm", "converted.cm", "converted.nopad.cm",
                "seed17a.filt", "seed17b.filt",
                "cmbuild.log") {
    if (-e "$tmppfx.$suffix") { unlink "$tmppfx.$suffix"; }
}
exit 0;

######################################################################
# Helper subs
######################################################################

sub first_line
{
    my ($file) = @_;
    open(FH, "<", $file) || die "FAIL: unable to open $file\n";
    my $line = <FH>;
    close FH;
    chomp $line if defined $line;
    return (defined $line) ? $line : "";
}

# True if file contains a header line "<tag> <value>" where <value>
# is the first whitespace-separated token after <tag>. Matches e.g.
# "P7NODEPAD yes" or "P7NODEPAD no".
sub file_has_header
{
    my ($file, $tag, $value) = @_;
    open(FH, "<", $file) || die "FAIL: unable to open $file\n";
    while (my $line = <FH>) {
        chomp $line;
        if ($line =~ /^\Q$tag\E\s+(\S+)/) {
            close FH;
            return ($1 eq $value) ? 1 : 0;
        }
    }
    close FH;
    return 0;
}

# Return the first node line starting with "[ MATL", trimmed leading
# whitespace preserved so the caller can do a whitespace split.
sub first_matl_line
{
    my ($file) = @_;
    open(FH, "<", $file) || die "FAIL: unable to open $file\n";
    while (my $line = <FH>) {
        chomp $line;
        if ($line =~ /\[\s*MATL\b/) { close FH; return $line; }
    }
    close FH;
    return undef;
}

# Strip lines that embed timestamps or command-line info (COM, DATE
# in both the CM header and the filter-HMM section) so two independent
# runs can be byte-compared.
sub strip_com_lines
{
    my ($in, $out) = @_;
    open(IN,  "<", $in)  || die "FAIL: unable to open $in\n";
    open(OUT, ">", $out) || die "FAIL: unable to open $out for write\n";
    while (my $line = <IN>) {
        next if $line =~ /^COM\s/;
        next if $line =~ /^DATE\s/;    # CM header timestamp
        next if $line =~ /^DATE\b/;    # filter-HMM section timestamp ("DATE  ...")
        print OUT $line;
    }
    close IN;
    close OUT;
}

# Pull the trailing 2 fields from every MAT* node line. Used to
# byte-compare pad columns independent of everything else.
sub extract_pad_columns
{
    my ($file) = @_;
    my $out = "";
    open(FH, "<", $file) || die "FAIL: unable to open $file\n";
    while (my $line = <FH>) {
        chomp $line;
        next unless $line =~ /\[\s*MAT[PLR]\b/;
        my @f = split /\s+/, $line;
        shift @f if $f[0] eq "";
        next if scalar(@f) < 12;
        # node-type indicator (e.g. MATL) is at idx 1, node index at
        # idx 2, so we prefix it into the extracted line so mismatches
        # are easy to read:
        $out .= "$f[1] $f[2] $f[10] $f[11]\n";
    }
    close FH;
    return $out;
}
