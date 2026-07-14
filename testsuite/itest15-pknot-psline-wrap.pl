#! /usr/bin/perl

# Regression test for the cmsearch pknot PS line (pair-status: '=' maintained,
# '$' covarying, 'x' broken) under DISPLAY LINE WRAPPING.
#
# The single riskiest property of the PS line is that the =/$/x marks stay
# column-aligned under the CS pknot letters when the cmsearch alignment WRAPS
# across multiple display blocks. Every other pknot test (test-010, itest13)
# uses clen 34-35 models that fit in one block, and cmsearch's --textw floor is
# 120, so a small model CANNOT exercise wrapping. This test builds a wide
# (clen 129) model whose pknot pair is long-range -- the open stem (AAAA, near
# display col 6) and the close stem (aaaa, near display col 126) land in
# DIFFERENT wrapped blocks at --textw 120 -- and asserts that in EACH block the
# PS-line marks occur at EXACTLY the column offsets of the pknot letters on that
# block's CS line.
#
# The assertion is non-vacuous: it compares per-column offset SETS (PS non-blank
# marks vs CS [A-Za-z] pknot letters) for byte-exact equality, so a one-column
# shift, a dropped mark, or a mark in the wrong block all FAIL. It also requires
# the alignment to actually wrap (>=2 blocks) and the open/close ends to land in
# different blocks, so the test cannot silently degrade to the already-covered
# single-block case.
#
# Targets:
#   maintained -- consensus sequence; pknot G...C pair stays WC -> marks are '='
#   broken     -- pknot close mutated C->G (G...G non-canonical) -> marks are 'x'
# ('$' covary is a separately-noted narrow case and is not exercised here.)
#
# Usage:    ./itest15-pknot-psline-wrap.pl <builddir> <srcdir> <tmpfile prefix>
# Example:  ./itest15-pknot-psline-wrap.pl ..         ..       tmpfoo

BEGIN {
    $builddir = shift;
    $srcdir   = shift;
    $tmppfx   = shift;
}

use strict;
use warnings;

our ($builddir, $srcdir, $tmppfx);

# Verify executables and the checked-in wide pknot seed.
if (! -x "$builddir/src/cmbuild")              { die "FAIL: didn't find cmbuild binary in $builddir/src\n";     }
if (! -x "$builddir/src/cmcalibrate")          { die "FAIL: didn't find cmcalibrate binary in $builddir/src\n"; }
if (! -x "$builddir/src/cmsearch")             { die "FAIL: didn't find cmsearch binary in $builddir/src\n";    }
if (! -e "$srcdir/testsuite/wide-pknot.sto")   { die "FAIL: didn't find wide-pknot.sto in $srcdir/testsuite\n"; }

my $seed = "$srcdir/testsuite/wide-pknot.sto";

my @tmpsuffixes = ("cm", "fa", "search.out", "log");
sub cleanup { for my $s (@tmpsuffixes) { unlink "$tmppfx.$s" if -e "$tmppfx.$s"; } }
cleanup();

######################################################################
# Step 1: build the wide pknot CM and confirm it is a pknot model that
# is wide enough to force wrapping at --textw 120 (clen >= 130 margin is
# not required; clen 129 already wraps -- assert clen >= 121 so the
# close stem is guaranteed past the first 120-col block, and PKNOT yes).
######################################################################
run("$builddir/src/cmbuild --wnone -F $tmppfx.cm $seed > $tmppfx.log 2>&1",
    "step 1: cmbuild failed");
if (! file_has_header("$tmppfx.cm", "PKNOT", "yes")) { die "FAIL: step 1: built CM missing 'PKNOT yes'\n"; }
my $clen = header_value("$tmppfx.cm", "CLEN");
if (!defined $clen || $clen < 121) { die "FAIL: step 1: CLEN is '" . (defined $clen ? $clen : "undef") . "', need >= 121 to force wrapping at --textw 120\n"; }

######################################################################
# Step 2: calibrate (cmsearch refuses to run without E-value params).
# Tiny -L keeps it fast; a fixed --seed keeps it reproducible. The
# calibration RNG only affects E-values, not the CYK alignment that
# drives the PS/CS columns, so this does not affect the assertion.
######################################################################
run("$builddir/src/cmcalibrate -L 0.05 --seed 19 $tmppfx.cm > $tmppfx.log 2>&1",
    "step 2: cmcalibrate failed");

######################################################################
# Step 3: write the two engineered search targets.
#   maintained : consensus -> pknot pair G...C is WC -> '=' marks
#   broken     : pknot close CCCC -> GGGG -> G...G non-canonical -> 'x' marks
# The nested stem (cols 1-5 / 121-125) is left intact in both so the hit
# still aligns full-length.
######################################################################
my $maint  = ("G" x 9) . ("A" x 111) . ("C" x 9);              # GGGGG GGGG ...A... CCCCC CCCC
my $broken = ("G" x 9) . ("A" x 111) . ("C" x 5) . ("G" x 4);  # ...                CCCCC GGGG  <- pknot close broken
open(my $fafh, ">", "$tmppfx.fa") || die "FAIL: step 3: unable to write $tmppfx.fa\n";
print $fafh ">maintained\n$maint\n>broken\n$broken\n";
close $fafh;

######################################################################
# Step 4: search at --textw 120 (forces wrapping) with --max -T 5 (so the
# engineered targets hit) and --cpu 0 (deterministic, single-threaded).
######################################################################
run("$builddir/src/cmsearch --max -T 5 --textw 120 --cpu 0 $tmppfx.cm $tmppfx.fa > $tmppfx.search.out 2>&1",
    "step 4: cmsearch failed");

######################################################################
# Step 5: per-target assertions on the wrapped alignment blocks.
######################################################################
my %blocks = parse_psline_blocks("$tmppfx.search.out");   # target -> [ {ps=>str, cs=>str}, ... ]

check_target(\%blocks, "maintained", "=");
check_target(\%blocks, "broken",     "x");

######################################################################
print "ok\n";
cleanup();
exit 0;

######################################################################
# Assertion driver for one target.
######################################################################
sub check_target {
    my ($blocks, $target, $expect_mark) = @_;
    my $blks = $blocks->{$target};
    if (!defined $blks) { die "FAIL: no '>> $target' alignment found in cmsearch output\n"; }

    # (a) Must actually wrap: >= 2 display blocks. Otherwise the test would
    # silently degrade to the already-covered single-block case.
    if (scalar(@$blks) < 2) {
        die "FAIL: target '$target': alignment did not wrap (got " . scalar(@$blks) . " block(s), need >= 2). --textw 120 should wrap a clen-129 model.\n";
    }

    my $saw_upper_block = -1;   # index of block holding the open (uppercase) pknot letters
    my $saw_lower_block = -1;   # index of block holding the close (lowercase) pknot letters
    my $total_marks = 0;

    for (my $b = 0; $b < scalar(@$blks); $b++) {
        my $ps = $blks->[$b]{ps};
        my $cs = $blks->[$b]{cs};

        # (b) THE crux: per-column offset sets must match exactly.
        my @cs_pk = offsets_matching($cs, qr/[A-Za-z]/);  # pknot letters on CS
        my @ps_mk = offsets_matching($ps, qr/\S/);         # non-blank marks on PS
        my $cs_key = join(",", @cs_pk);
        my $ps_key = join(",", @ps_mk);
        if ($cs_key ne $ps_key) {
            die "FAIL: target '$target' block " . ($b+1) . ": PS marks are NOT column-aligned with CS pknot letters.\n"
              . "  CS pknot-letter offsets: [$cs_key]\n"
              . "  PS mark         offsets: [$ps_key]\n"
              . "  CS: |$cs|\n  PS: |$ps|\n";
        }

        # (c) Every mark in this block must be the expected glyph.
        for my $i (@ps_mk) {
            my $m = substr($ps, $i, 1);
            if ($m ne $expect_mark) {
                die "FAIL: target '$target' block " . ($b+1) . ": mark at offset $i is '$m', expected '$expect_mark'\n";
            }
            $total_marks++;
        }

        # Track which block carries the open (uppercase) vs close (lowercase) end.
        if (grep { substr($cs, $_, 1) =~ /[A-Z]/ } @cs_pk) { $saw_upper_block = $b if $saw_upper_block < 0; }
        if (grep { substr($cs, $_, 1) =~ /[a-z]/ } @cs_pk) { $saw_lower_block = $b if $saw_lower_block < 0; }
    }

    # (d) There must actually BE marks (guards against a vacuous empty match).
    if ($total_marks == 0) { die "FAIL: target '$target': no PS marks found at all\n"; }

    # (e) Open and close ends must land in DIFFERENT blocks -> the pknot pair
    # genuinely spans the wrap boundary (the property under test). This also
    # rules out "all marks crammed into one block".
    if ($saw_upper_block < 0) { die "FAIL: target '$target': no open-stem (uppercase) pknot letters seen on any CS line\n"; }
    if ($saw_lower_block < 0) { die "FAIL: target '$target': no close-stem (lowercase) pknot letters seen on any CS line\n"; }
    if ($saw_upper_block == $saw_lower_block) {
        die "FAIL: target '$target': open and close pknot ends are in the SAME block ($saw_upper_block); the pair did not span a wrap boundary, so wrapping is untested\n";
    }
}

######################################################################
# Parse cmsearch human-readable output into per-target lists of
# {ps,cs} display blocks. The PS line ends in ' PS', the CS line is the
# immediately following ' CS' line; both share the same field width
# before the trailer, so column offsets are directly comparable once the
# trailing ' PS'/' CS' is stripped (and only that -- internal spacing is
# preserved). Within each target, blocks are returned in display order.
######################################################################
sub parse_psline_blocks {
    my ($file) = @_;
    my %out;
    open(my $fh, "<", $file) || die "FAIL: unable to open $file\n";
    my @lines = <$fh>;
    close $fh;
    chomp @lines;

    my $cur;                    # current target name
    for (my $i = 0; $i < scalar(@lines); $i++) {
        my $line = $lines[$i];
        if ($line =~ /^>>\s+(\S+)/) { $cur = $1; $out{$cur} = []; next; }
        next unless defined $cur;
        # A PS line (ends with the ' PS' trailer) paired with the next ' CS' line.
        if ($line =~ / PS$/ && $i + 1 < scalar(@lines) && $lines[$i+1] =~ / CS$/) {
            (my $ps = $line)        =~ s/ PS$//;
            (my $cs = $lines[$i+1]) =~ s/ CS$//;
            push @{$out{$cur}}, { ps => $ps, cs => $cs };
            $i++;               # consume the CS line
        }
    }
    return %out;
}

# Sorted list of 0-based offsets where the string matches the given regex
# (applied per character).
sub offsets_matching {
    my ($str, $re) = @_;
    my @pos;
    for (my $i = 0; $i < length($str); $i++) {
        push @pos, $i if substr($str, $i, 1) =~ $re;
    }
    return @pos;
}

######################################################################
# Generic helpers (mirrors itest13 style).
######################################################################
sub run {
    my ($cmd, $msg) = @_;
    my $out = `$cmd`;
    if ($? != 0) { die "FAIL: $msg\n$out\n"; }
}

# First whitespace token after a header tag, or undef.
sub header_value {
    my ($file, $tag) = @_;
    open(my $fh, "<", $file) || die "FAIL: unable to open $file\n";
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ /^\Q$tag\E\s+(\S+)/) { close $fh; return $1; }
    }
    close $fh;
    return undef;
}

# True if file has a header line "<tag> <value>" (value = first token after tag).
sub file_has_header {
    my ($file, $tag, $value) = @_;
    my $v = header_value($file, $tag);
    return (defined $v && $v eq $value) ? 1 : 0;
}
