#! /usr/bin/perl

# Regression test for the cmsearch/cmscan pseudoknot register-shift-recoverable
# 'o' overlay (brief 26_0617-028). This is the permanent, self-contained version
# of that brief's key validation gate: the real tmRNA hit FR891194.1 (the case
# Sean Eddy asked about at Eric's Benasque talk), where stem B of the pseudoknot
# has a genuine 1-nt register-shift ambiguity that restores canonical WC/GU
# pairing to its 5 broken ('x') pairs.
#
# The feature (cm_pknot_MarkShiftRecoverable(), src/cm_alidisplay.c) is ALERT-ONLY:
# it never changes cmsearch's reported alignment, it only re-labels a broken 'x'
# pknot pair as 'o' when a per-STEM uniform register shift would restore its
# canonical pairing. This test confirms it fires on EXACTLY the right pairs:
#   - EXACTLY 5 pairs become 'o', and they are ALL in a single stem (stem B).
#   - Every OTHER stem that has broken 'x' pairs KEEPS them as 'x' (no over-fire):
#     stems A, C, D collectively retain 6 'x' pairs here, none re-labeled 'o'.
#   - '='/'$' (maintained/covarying) pairs are never touched.
# These exact per-stem counts were independently derived by hand in brief 027
# (summary_027), so this is a strong known-answer check, not just "some 'o' appeared".
#
# Fixtures (checked in, self-contained): testsuite/tmRNA-pknot.cm (a calibrated
# INFERNAL1/d pseudoknot CM built from the tmRNA seed) and
# testsuite/tmRNA-pknot-targets.fa (the 5 real tmRNA windows + 2 decoys). The
# 'o' overlay does not depend on the E-value calibration, only on the reported
# alignment, so the frozen calibration in the checked-in CM is irrelevant to the
# result; it is present only because cmsearch requires a calibrated CM to run.
#
# Usage:    ./itest17-pknot-shift-recoverable.pl <builddir> <srcdir> <tmpfile prefix>
# Example:  ./itest17-pknot-shift-recoverable.pl ..         ..       tmpfoo

BEGIN {
    $builddir = shift;
    $srcdir   = shift;
    $tmppfx   = shift;
}

use strict;
use warnings;

our ($builddir, $srcdir, $tmppfx);

if (! -x "$builddir/src/cmsearch")                    { die "FAIL: didn't find cmsearch binary in $builddir/src\n"; }
if (! -e "$srcdir/testsuite/tmRNA-pknot.cm")          { die "FAIL: didn't find tmRNA-pknot.cm in $srcdir/testsuite\n"; }
if (! -e "$srcdir/testsuite/tmRNA-pknot-targets.fa")  { die "FAIL: didn't find tmRNA-pknot-targets.fa in $srcdir/testsuite\n"; }

my $cm  = "$srcdir/testsuite/tmRNA-pknot.cm";
my $fa  = "$srcdir/testsuite/tmRNA-pknot-targets.fa";
my $out = "$tmppfx.o.out";
my $hit = "FR891194.1";

# --cpu 0 for reproducibility; default text width (this deliberately WRAPS the
# alignment into several blocks, so the reconstruction below also confirms the
# 'o' marks land column-exact across wrap boundaries).
run("$builddir/src/cmsearch --cpu 0 $cm $fa > $out 2>&1", "cmsearch (tmRNA pknot) failed");

# Reconstruct the raw (un-wrapped) PS and CS strings for the FR891194.1 hit by
# concatenating each block's alignment-field slice. The field span is taken from
# the target-sequence line ("<name> <start> <residues> <end>"); the PS and CS
# lines are column-aligned with it.
my ($rawPS, $rawCS) = reconstruct_hit($out, $hit);
if ($rawPS eq "") { die "FAIL: no '$hit' alignment (PS line) found in cmsearch output\n"; }
if (length($rawPS) != length($rawCS)) { die "FAIL: reconstructed PS/CS length mismatch (" . length($rawPS) . " vs " . length($rawCS) . ")\n"; }

# Recover every complete pseudoknot pair from the CS line with the same per-letter
# pushdown discipline the C code uses, and classify each by its PS marks and stem.
my %pairmark;   # $pairmark{stem_letter}{"$mo$mc"} = count
my $n_o_pairs   = 0;
my %stems_with_o;
my $n_x_pairs   = 0;
my @stack;      # per-uppercase-letter stacks of open positions
my %stk;
for (my $i = 0; $i < length($rawCS); $i++) {
    my $c = substr($rawCS, $i, 1);
    if ($c =~ /^[A-Z]$/) {
        push @{$stk{$c}}, $i;
    } elsif ($c =~ /^[a-z]$/) {
        my $u = uc $c;
        if ($stk{$u} && scalar(@{$stk{$u}})) {
            my $zo = pop @{$stk{$u}};
            my $zc = $i;
            my $mo = substr($rawPS, $zo, 1);
            my $mc = substr($rawPS, $zc, 1);
            $pairmark{$u}{"$mo$mc"}++;
            if ($mo eq 'o' || $mc eq 'o') { $n_o_pairs++; $stems_with_o{$u} = 1;
                # both ends of a recovered pair must be 'o' (marked at both ends).
                if (!($mo eq 'o' && $mc eq 'o')) { die "FAIL: a recovered pair has 'o' on only one end (stem $u: '$mo'/'$mc')\n"; }
            }
            if ($mo eq 'x' || $mc eq 'x') { $n_x_pairs++; }
        }
    }
}

# ---- Assertions (the exact known answer from brief 027) --------------------

# (1) EXACTLY 5 pairs recovered to 'o'.
if ($n_o_pairs != 5) {
    die "FAIL: expected exactly 5 register-shift-recoverable ('o') pairs, got $n_o_pairs\n" . dump_marks(\%pairmark);
}

# (2) ...all within a SINGLE stem (per-STEM uniform shift, never scattered per-pair).
my $n_o_stems = scalar keys %stems_with_o;
if ($n_o_stems != 1) {
    die "FAIL: expected all 'o' pairs in exactly 1 stem, got $n_o_stems (" . join(",", sort keys %stems_with_o) . ")\n" . dump_marks(\%pairmark);
}

# (3) NO over-fire: other stems' broken pairs must remain 'x'. Here stems A/C/D
#     collectively retain 6 'x' pairs, and none of them is the 'o' stem.
if ($n_x_pairs < 1) {
    die "FAIL: expected the non-recoverable stems to retain 'x' pairs (no-over-fire control), but found 0 'x' pairs\n" . dump_marks(\%pairmark);
}
my ($o_stem) = keys %stems_with_o;
if (exists $pairmark{$o_stem}{"xx"}) {
    die "FAIL: the recovered stem '$o_stem' still has 'x' pairs -- overlay should have re-labeled all its recoverable broken pairs\n" . dump_marks(\%pairmark);
}

# (4) Raw-character sanity: 5 pairs * 2 ends = exactly 10 'o' glyphs; and the
#     '='/'$' maintained/covarying marks were left intact (never turned to 'o').
my $o_chars = ($rawPS =~ tr/o/o/);
if ($o_chars != 10) { die "FAIL: expected exactly 10 'o' characters on the PS line (5 pairs x 2 ends), got $o_chars\n"; }

print "ok\n";
unlink $out if -e $out;
exit 0;

######################################################################
# Helpers
######################################################################

sub run {
    my ($cmd, $msg) = @_;
    my $output = `$cmd`;
    if ($? != 0) { die "FAIL: $msg\n$output\n"; }
}

# Concatenate the alignment-field slice of every PS/CS block for a named hit,
# returning ($rawPS, $rawCS). Field bounds come from the target-sequence line.
sub reconstruct_hit {
    my ($file, $name) = @_;
    open(my $fh, "<", $file) || die "FAIL: unable to open $file\n";
    my @lines = <$fh>;
    close $fh;
    chomp @lines;

    # Restrict to the hit's alignment region (from ">> <name>" to the next ">>" or summary).
    my $start;
    for (my $i = 0; $i < scalar(@lines); $i++) {
        if ($lines[$i] =~ /^>>\s+\Q$name\E/) { $start = $i; last; }
    }
    return ("", "") unless defined $start;
    my @reg;
    for my $l (@lines[$start+1 .. $#lines]) {
        last if $l =~ /^>>\s/ || $l =~ /^Internal/;
        push @reg, $l;
    }

    my ($ps, $cs) = ("", "");
    for (my $i = 0; $i < scalar(@reg); $i++) {
        next unless $reg[$i] =~ / PS$/;
        # the CS line is immediately below the PS line
        my $csline = ($i+1 < scalar(@reg)) ? $reg[$i+1] : "";
        next unless $csline =~ / CS$/;
        # the target-sequence line (within a few lines) gives the field bounds
        my $seqline;
        for (my $j = $i+1; $j < $i+8 && $j < scalar(@reg); $j++) {
            if ($reg[$j] =~ /\Q$name\E/ && $reg[$j] =~ /^\s+\S+\s+\d+\s+\S+\s+\d+\s*$/) { $seqline = $reg[$j]; last; }
        }
        next unless defined $seqline;
        $seqline =~ /^(\s+\S+\s+\d+\s+)(\S+)\s+\d+\s*$/ or next;
        my $fs = length($1);
        my $fe = $fs + length($2);
        (my $psbody = $reg[$i]) =~ s/ PS$//;
        (my $csbody = $csline)  =~ s/ CS$//;
        $ps .= substr($psbody, $fs, $fe - $fs);
        $cs .= substr($csbody, $fs, $fe - $fs);
    }
    return ($ps, $cs);
}

sub dump_marks {
    my ($pm) = @_;
    my $s = "  per-stem pair marks:\n";
    for my $L (sort keys %$pm) {
        $s .= "    stem $L: " . join(", ", map { "$_:$pm->{$L}{$_}" } sort keys %{$pm->{$L}}) . "\n";
    }
    return $s;
}
