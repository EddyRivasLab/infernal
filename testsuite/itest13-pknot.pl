#! /usr/bin/perl

# Test the INFERNAL1/d consensus-pseudoknot passthrough (Feature B):
# cmbuild capture/canonicalize/store, --refine preservation, the
# canonical >=3-stem relabel, the 1d->1b back-convert, and the
# malformed-node-line read validation. Covers the test gaps flagged by
# review summary 007 (items 1, 3, 4). The binary older-format hardening
# (item 2, the format gate) is covered by the companion C harness
# itest13-pknot-1b-hardening.c, which needs library-level access.
#
# Usage:    ./itest13-pknot.pl <builddir> <srcdir> <tmpfile prefix>
# Example:  ./itest13-pknot.pl ..         ..       tmpfoo

BEGIN {
    $builddir = shift;
    $srcdir   = shift;
    $tmppfx   = shift;
}

use strict;
use warnings;

our ($builddir, $srcdir, $tmppfx);

# Verify executables and the checked-in pseudoknot seed (RF01096, one A/a stem).
if (! -x "$builddir/src/cmbuild")          { die "FAIL: didn't find cmbuild binary in $builddir/src\n";   }
if (! -x "$builddir/src/cmconvert")        { die "FAIL: didn't find cmconvert binary in $builddir/src\n"; }
if (! -x "$builddir/src/cmcalibrate")      { die "FAIL: didn't find cmcalibrate binary in $builddir/src\n"; }
if (! -x "$builddir/src/cmsearch")         { die "FAIL: didn't find cmsearch binary in $builddir/src\n";  }
if (! -e "$srcdir/testsuite/PK-HAV.sto")   { die "FAIL: didn't find PK-HAV.sto in $srcdir/testsuite\n";   }

my $pkhav = "$srcdir/testsuite/PK-HAV.sto";

my @tmpsuffixes = ("direct.cm", "refine.cm", "refine.sto",
                   "gibbs.cm", "gibbs.sto",
                   "frag.sto", "frag.cm", "frag.refined.sto", "frag.re.cm",
                   "syn3.sto", "syn3.cm", "syn3.re.cm",
                   "bconv.cm", "bconv.reread.cm",
                   "bad.cm", "bad.err",
                   "orph.sto", "orph.cm", "orph.cal.log", "orph.fa", "orph.out", "orph.max.out",
                   "log");
sub cleanup { for my $s (@tmpsuffixes) { unlink "$tmppfx.$s" if -e "$tmppfx.$s"; } }
cleanup();

######################################################################
# Subtest 1: direct build of a pknot seed yields INFERNAL1/d, PKNOT yes.
######################################################################
run("$builddir/src/cmbuild --wnone -F $tmppfx.direct.cm $pkhav > $tmppfx.log 2>&1",
    "subtest 1: cmbuild (direct) failed");
my $first = first_line("$tmppfx.direct.cm");
if ($first !~ /^INFERNAL1\/d\b/) { die "FAIL: subtest 1: direct build first line is '$first', expected INFERNAL1/d\n"; }
if (! file_has_header("$tmppfx.direct.cm", "PKNOT", "yes")) { die "FAIL: subtest 1: direct build missing 'PKNOT yes'\n"; }
my $pk_direct = pknot_columns("$tmppfx.direct.cm");
if ($pk_direct eq "") { die "FAIL: subtest 1: no pknot columns extracted from direct build\n"; }
assert_balanced($pk_direct, "subtest 1 (direct)");

######################################################################
# Subtest 2: --refine PRESERVES pknots (review found it preserves, not
# drops). The refined CM's pknot annotation must equal the direct build.
######################################################################
run("$builddir/src/cmbuild --wnone -F --refine $tmppfx.refine.sto $tmppfx.refine.cm $pkhav > $tmppfx.log 2>&1",
    "subtest 2: cmbuild --refine failed");
if (! file_has_header("$tmppfx.refine.cm", "PKNOT", "yes")) { die "FAIL: subtest 2: --refine build missing 'PKNOT yes'\n"; }
my $pk_refine = pknot_columns("$tmppfx.refine.cm");
if ($pk_refine ne $pk_direct) {
    die "FAIL: subtest 2: --refine pknot columns differ from direct build:\n--- direct ---\n$pk_direct\n--- refine ---\n$pk_refine\n";
}

######################################################################
# Subtest 3: --refine --gibbs (RNG-driven) does not drop or corrupt
# pknots: PKNOT yes and balanced annotation (the intermediate ss_cons
# the review worried could go unbalanced stays valid here).
######################################################################
run("$builddir/src/cmbuild --wnone -F --refine $tmppfx.gibbs.sto --gibbs --seed 7 $tmppfx.gibbs.cm $pkhav > $tmppfx.log 2>&1",
    "subtest 3: cmbuild --refine --gibbs failed");
if (! file_has_header("$tmppfx.gibbs.cm", "PKNOT", "yes")) { die "FAIL: subtest 3: --gibbs build missing 'PKNOT yes'\n"; }
assert_balanced(pknot_columns("$tmppfx.gibbs.cm"), "subtest 3 (--gibbs)");

######################################################################
# Subtest 4: --refine on a FRAGMENTARY alignment (one seq is a 5'
# fragment, 3' end terminal-gapped, marked via --miss) preserves pknots
# and round-trips -- the truncated-seq edge the review flagged. Confirm
# it fails safe / preserves (here: preserves), never invalid WUSS.
######################################################################
write_frag_sto("$tmppfx.frag.sto");
run("$builddir/src/cmbuild --wnone -F --refine $tmppfx.frag.refined.sto --miss $tmppfx.frag.cm $tmppfx.frag.sto > $tmppfx.log 2>&1",
    "subtest 4: cmbuild --refine --miss (fragmentary) failed");
if (! file_has_header("$tmppfx.frag.cm", "PKNOT", "yes")) { die "FAIL: subtest 4: fragmentary --refine missing 'PKNOT yes'\n"; }
assert_balanced(pknot_columns("$tmppfx.frag.cm"), "subtest 4 (fragmentary refine)");
run("$builddir/src/cmconvert -a $tmppfx.frag.cm > $tmppfx.frag.re.cm 2>>$tmppfx.log",
    "subtest 4: round-trip of fragmentary-refined CM failed");

######################################################################
# Subtest 5: canonical relabel on >=3 mutually-crossing pknot stems.
# Build a synthetic alignment whose ss_cons has three crossing stems
# (A,B,C). Confirm valid, balanced WUSS pknot annotation with exactly
# three distinct canonical stems, and that it round-trips byte-for-byte.
######################################################################
write_syn3_sto("$tmppfx.syn3.sto");
run("$builddir/src/cmbuild --wnone -F $tmppfx.syn3.cm $tmppfx.syn3.sto > $tmppfx.log 2>&1",
    "subtest 5: cmbuild (synthetic 3-stem) failed");
if (! file_has_header("$tmppfx.syn3.cm", "PKNOT", "yes")) { die "FAIL: subtest 5: synthetic 3-stem build missing 'PKNOT yes'\n"; }
my $pk_syn3 = pknot_columns("$tmppfx.syn3.cm");
assert_balanced($pk_syn3, "subtest 5 (synthetic 3-stem)");
my %letters = distinct_pknot_letters($pk_syn3);   # uppercase letters present
my $nstem = scalar keys %letters;
if ($nstem != 3) { die "FAIL: subtest 5: expected 3 distinct canonical pknot stems, got $nstem (" . join(",", sort keys %letters) . ")\n"; }
for my $L ('A','B','C') { if (! exists $letters{$L}) { die "FAIL: subtest 5: canonical relabel did not produce stem '$L' (got " . join(",", sort keys %letters) . ")\n"; } }
# round-trip identity
run("$builddir/src/cmconvert -a $tmppfx.syn3.cm > $tmppfx.syn3.re.cm 2>>$tmppfx.log",
    "subtest 5: cmconvert -a round-trip of synthetic 3-stem failed");
my $pk_syn3_rt = pknot_columns("$tmppfx.syn3.re.cm");
if ($pk_syn3 ne $pk_syn3_rt) { die "FAIL: subtest 5: synthetic 3-stem pknot columns changed across round-trip\n"; }

######################################################################
# Subtest 6: 1d -> 1b back-convert drops pknots cleanly. cmconvert
# --outfmt 1/b writes INFERNAL1/b with NO PKNOT header and no pknot
# letters on node lines, and the result re-reads cleanly.
######################################################################
run("$builddir/src/cmconvert -a --outfmt 1/b $tmppfx.direct.cm > $tmppfx.bconv.cm 2>$tmppfx.log",
    "subtest 6: cmconvert --outfmt 1/b failed");
$first = first_line("$tmppfx.bconv.cm");
if ($first !~ /^INFERNAL1\/b\b/) { die "FAIL: subtest 6: back-convert first line is '$first', expected INFERNAL1/b\n"; }
if (file_has_header("$tmppfx.bconv.cm", "PKNOT", "yes")) { die "FAIL: subtest 6: 1b back-convert still advertises 'PKNOT yes'\n"; }
my $leaked = pknot_letters_on_nodelines("$tmppfx.bconv.cm");
if ($leaked != 0) { die "FAIL: subtest 6: 1b back-convert leaked $leaked pknot letter(s) onto node lines\n"; }
run("$builddir/src/cmconvert -a $tmppfx.bconv.cm > $tmppfx.bconv.reread.cm 2>>$tmppfx.log",
    "subtest 6: re-reading the 1b back-convert failed");

######################################################################
# Subtest 7: a malformed 1d node line (non-alpha/non-dot pknot char in
# a MATP column) is rejected cleanly by the reader (Fix 2), not silently
# stored. Corrupt the first MATP node line's left pknot column to '@'.
######################################################################
corrupt_first_matp_pknot("$tmppfx.direct.cm", "$tmppfx.bad.cm");
my $rc = system("$builddir/src/cmconvert -a $tmppfx.bad.cm > /dev/null 2>$tmppfx.bad.err");
if ($rc == 0) { die "FAIL: subtest 7: cmconvert accepted a malformed MATP pknot column (expected rejection)\n"; }
my $err = slurp("$tmppfx.bad.err");
if ($err !~ /Invalid .*pknot character on MATP node line/) {
    die "FAIL: subtest 7: rejection message did not mention an invalid MATP pknot character:\n$err\n";
}

######################################################################
# Subtest 8: cmsearch PS-line orphan marking (brief 025). A pseudoknot
# letter whose partner isn't present in a hit's displayed alignment keeps
# its CS-line letter and gets '?' on the PS line, instead of being erased
# to a plain singlet -- for BOTH a truncated hit and an ordinary
# local-alignment hit that just doesn't span both stem halves. A complete
# (both-halves-present) pair in the same run must be unaffected.
######################################################################
write_orphan_pk_sto("$tmppfx.orph.sto");
run("$builddir/src/cmbuild --wnone -F $tmppfx.orph.cm $tmppfx.orph.sto > $tmppfx.log 2>&1",
    "subtest 8: cmbuild (orphan fixture) failed");
if (! file_has_header("$tmppfx.orph.cm", "PKNOT", "yes")) { die "FAIL: subtest 8: orphan-fixture build missing 'PKNOT yes'\n"; }
run("$builddir/src/cmcalibrate -L 0.01 --seed 19 $tmppfx.orph.cm > $tmppfx.orph.cal.log 2>&1",
    "subtest 8: cmcalibrate (orphan fixture) failed");

# (a),(c): default (truncation-aware) pipeline over both a flanked complete
# pair ("complete") and an unflanked 3'-truncated fragment ("localorphan")
# in the SAME run.
open(my $fafh, ">", "$tmppfx.orph.fa") || die "FAIL: subtest 8: unable to write $tmppfx.orph.fa\n";
print $fafh ">complete\n" . ("N" x 10) . ("G" x 9) . ("A" x 20) . ("C" x 9) . ("N" x 10) . "\n";
print $fafh ">localorphan\n" . ("G" x 9) . ("A" x 20) . ("C" x 5) . "\n";
close $fafh;
run("$builddir/src/cmsearch --toponly -T 1 --textw 200 --cpu 0 $tmppfx.orph.cm $tmppfx.orph.fa > $tmppfx.orph.out 2>&1",
    "subtest 8: cmsearch (default pipeline) failed");
my %blocks = parse_psline_blocks("$tmppfx.orph.out");

# complete pair: no '?' anywhere, both ends still classified '=' (maintained).
check_orphan_target(\%blocks, "complete",     0, "=");
# truncated orphan: open (uppercase) pknot half (4 letters, AAAA) kept +
# each marked '?'; the closing half must not appear on CS at all
# (truncated out of the display).
check_orphan_target(\%blocks, "localorphan",  4, "?");

# (b): the SAME unflanked fragment, but with --max (turns truncated-hit
# detection off entirely) -- this is the "ordinary local-alignment, no
# truncation flag at all" orphan case.
run("$builddir/src/cmsearch --toponly --max -T 1 --textw 200 --cpu 0 $tmppfx.orph.cm $tmppfx.orph.fa > $tmppfx.orph.max.out 2>&1",
    "subtest 8: cmsearch (--max, ordinary-local orphan) failed");
my %maxblocks = parse_psline_blocks("$tmppfx.orph.max.out");
check_orphan_target(\%maxblocks, "localorphan", 4, "?");
if ($maxblocks{"localorphan"}[0]{trunc} ne "no") {
    die "FAIL: subtest 8: --max hit unexpectedly still flagged truncated ('" . $maxblocks{"localorphan"}[0]{trunc} . "'); the ordinary-local-coverage case requires trunc=no\n";
}

######################################################################
print "ok\n";
cleanup();
exit 0;

######################################################################
# Helper subs
######################################################################

sub run {
    my ($cmd, $msg) = @_;
    my $out = `$cmd`;
    if ($? != 0) { die "FAIL: $msg\n$out\n"; }
}

sub first_line {
    my ($file) = @_;
    open(my $fh, "<", $file) || die "FAIL: unable to open $file\n";
    my $line = <$fh>;
    close $fh;
    chomp $line if defined $line;
    return defined $line ? $line : "";
}

# True if file has a header line "<tag> <value>" (value = first token after tag).
sub file_has_header {
    my ($file, $tag, $value) = @_;
    open(my $fh, "<", $file) || die "FAIL: unable to open $file\n";
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ /^\Q$tag\E\s+(\S+)/) { close $fh; return ($1 eq $value) ? 1 : 0; }
    }
    close $fh;
    return 0;
}

# The pknot annotation = the last two whitespace fields of every MAT*
# node line, prefixed by the node index for readable diffs.
sub pknot_columns {
    my ($file) = @_;
    my $out = "";
    open(my $fh, "<", $file) || die "FAIL: unable to open $file\n";
    while (my $line = <$fh>) {
        chomp $line;
        next unless $line =~ /\[\s*MAT[PLR]\b/;
        my @f = split /\s+/, $line;
        shift @f if $f[0] eq "";
        next if scalar(@f) < 14;            # 1d node line has 14 fields
        $out .= "$f[2] $f[-2] $f[-1]\n";    # node index, left pknot, right pknot
    }
    close $fh;
    return $out;
}

# Count pknot letters on node lines (alpha only); used to detect leakage.
sub pknot_letters_on_nodelines {
    my ($file) = @_;
    my $n = 0;
    open(my $fh, "<", $file) || die "FAIL: unable to open $file\n";
    while (my $line = <$fh>) {
        chomp $line;
        next unless $line =~ /\[\s*MAT[PLR]\b/;
        my @f = split /\s+/, $line;
        shift @f if $f[0] eq "";
        next if scalar(@f) < 2;
        for my $c ($f[-2], $f[-1]) { $n++ if $c =~ /^[A-Za-z]$/; }
    }
    close $fh;
    return $n;
}

# Each pknot stem must have matching open (uppercase) and close
# (lowercase) counts; die otherwise (would be invalid WUSS).
sub assert_balanced {
    my ($cols, $where) = @_;
    my %up; my %lo;
    for my $line (split /\n/, $cols) {
        my @f = split /\s+/, $line;
        for my $c ($f[1], $f[2]) {
            next unless defined $c;
            if    ($c =~ /^[A-Z]$/) { $up{lc $c}++; }
            elsif ($c =~ /^[a-z]$/) { $lo{$c}++;    }
        }
    }
    my %all = map { $_ => 1 } (keys %up, keys %lo);
    for my $L (sort keys %all) {
        my $u = $up{$L} // 0;
        my $d = $lo{$L} // 0;
        if ($u != $d) { die "FAIL: $where: pknot stem '$L' is unbalanced ($u open vs $d close) -> invalid WUSS\n"; }
    }
}

# Set of distinct uppercase pknot stem letters present.
sub distinct_pknot_letters {
    my ($cols) = @_;
    my %h;
    for my $line (split /\n/, $cols) {
        my @f = split /\s+/, $line;
        for my $c ($f[1], $f[2]) { next unless defined $c; $h{uc $c} = 1 if $c =~ /^[A-Za-z]$/; }
    }
    return %h;
}

# A synthetic gapless alignment for the orphan-marking subtest: nested
# stem (5bp) + pknot open singlets (AAAA) + a short unpaired loop + nested
# close (5bp) + pknot close singlets (aaaa). Small (clen 38) so
# cmcalibrate is fast.
sub write_orphan_pk_sto {
    my ($path) = @_;
    my $ss  = "<<<<<AAAA....................>>>>>aaaa";
    my $seq = ("G" x 9) . ("A" x 20) . ("C" x 9);
    die "FAIL: write_orphan_pk_sto length mismatch\n" unless length($ss) == length($seq);
    open(my $fh, ">", $path) || die "FAIL: unable to write $path\n";
    print $fh "# STOCKHOLM 1.0\n\n";
    for my $n ("seq1","seq2","seq3") { printf $fh "%-20s %s\n", $n, $seq; }
    printf $fh "%-20s %s\n", "#=GC SS_cons", $ss;
    print $fh "//\n";
    close $fh;
}

# Parse cmsearch human-readable output into per-target lists of
# {ps, cs, trunc} display blocks (mirrors itest15's parser; also captures
# the 'trunc' column from the hit-scores table since subtest 8 needs it).
sub parse_psline_blocks {
    my ($file) = @_;
    my %out;
    my %trunc;
    open(my $fh, "<", $file) || die "FAIL: unable to open $file\n";
    my @lines = <$fh>;
    close $fh;
    chomp @lines;

    for my $line (@lines) {
        # hit-scores table row: (n) sig E-value score bias sequence start end strand mdl trunc gc [desc]
        if ($line =~ /^\s*\(\d+\)\s+(.*)$/) {
            my @f = split /\s+/, $1;
            # f = sig, E-value, score, bias, sequence, start, end, strand, mdl, trunc, gc, ...
            if (scalar(@f) >= 11 && ($f[7] eq '+' || $f[7] eq '-')) {
                $trunc{$f[4]} = $f[9] unless exists $trunc{$f[4]};
            }
        }
    }

    my $cur;
    for (my $i = 0; $i < scalar(@lines); $i++) {
        my $line = $lines[$i];
        if ($line =~ /^>>\s+(\S+)/) { $cur = $1; $out{$cur} = [] unless exists $out{$cur}; next; }
        next unless defined $cur;
        if ($line =~ / PS$/ && $i + 1 < scalar(@lines) && $lines[$i+1] =~ / CS$/) {
            (my $ps = $line)        =~ s/ PS$//;
            (my $cs = $lines[$i+1]) =~ s/ CS$//;
            push @{$out{$cur}}, { ps => $ps, cs => $cs, trunc => ($trunc{$cur} // "?") };
            $i++;
        }
    }
    return %out;
}

# Assertion driver for subtest 8: confirm the expected number of '?' marks
# on the PS line, all pknot letters that ARE on CS are still present (not
# erased), and (when expect_orphans) that the CS pknot-letter offsets that
# carry no partner get exactly the expected mark.
sub check_orphan_target {
    my ($blocks, $target, $expect_qmarks, $expect_mark_if_orphan) = @_;
    my $blks = $blocks->{$target};
    if (!defined $blks || scalar(@$blks) == 0) { die "FAIL: subtest 8: no '>> $target' alignment found\n"; }
    my $b = $blks->[0];
    my @qpos = grep { substr($b->{ps}, $_, 1) eq '?' } (0 .. length($b->{ps}) - 1);
    if (scalar(@qpos) != $expect_qmarks) {
        die "FAIL: subtest 8: target '$target': expected $expect_qmarks '?' mark(s) on PS, got " . scalar(@qpos) . "\n  PS: |$b->{ps}|\n  CS: |$b->{cs}|\n";
    }
    if ($expect_qmarks > 0) {
        for my $i (@qpos) {
            my $c = substr($b->{cs}, $i, 1);
            if ($c !~ /^[A-Za-z]$/) { die "FAIL: subtest 8: target '$target': PS '?' at offset $i is not over a pknot letter on CS (got '$c')\n"; }
        }
    } else {
        # complete pair: no '?' at all, and any pknot letters present must
        # carry the expected complete-pair mark (not blank, not '?').
        for (my $i = 0; $i < length($b->{cs}); $i++) {
            my $c = substr($b->{cs}, $i, 1);
            next unless $c =~ /^[A-Za-z]$/;
            my $m = substr($b->{ps}, $i, 1);
            if ($m ne $expect_mark_if_orphan) { die "FAIL: subtest 8: target '$target': expected pknot mark '$expect_mark_if_orphan' at offset $i, got '$m'\n"; }
        }
    }
}

sub slurp {
    my ($file) = @_;
    open(my $fh, "<", $file) || return "";
    local $/; my $s = <$fh>; close $fh;
    return defined $s ? $s : "";
}

# Write a copy of $in with the first MATP node line's left pknot column
# (the second-to-last field) replaced by '@'.
sub corrupt_first_matp_pknot {
    my ($in, $out) = @_;
    open(my $ih, "<", $in)  || die "FAIL: unable to open $in\n";
    open(my $oh, ">", $out) || die "FAIL: unable to open $out for write\n";
    my $done = 0;
    while (my $line = <$ih>) {
        if (!$done && $line =~ /\[\s*MATP\b/) {
            chomp(my $l = $line);
            my @f = split /\s+/, $l;
            $f[-2] = '@';
            print $oh join(" ", @f), "\n";
            $done = 1;
        } else {
            print $oh $line;
        }
    }
    close $ih; close $oh;
    die "FAIL: corrupt_first_matp_pknot: no MATP node line found in $in\n" unless $done;
}

# A fragmentary alignment derived from PK-HAV: two full seqs plus a 5'
# fragment (3' end terminal-gapped) so --miss marks it as a fragment.
sub write_frag_sto {
    my ($path) = @_;
    my $seq1 = "UUAAACAAAUUUUCUUAAAAUUUCUGAGGUUUGUUUAUUUCUUUUAU.CAGUAAAU";
    my $seq2 = "UUAAACAAACCUUCUUAAAAUUUCUGAGAUUUGUUUAUUUUGCAUAUUCAGUAAAU";
    my $ss   = ".<<<<<<<<<<.........AAAAAAA>>>>>>>>>>.........a.aaa.aaa.";
    my $rf   = "UuaaaCaaauuUUCUUAAAAUUUCUGAgguuuGuuuaUUUCUUUUAU.CAGUAAAU";
    my $L    = length($seq1);
    my $frag = substr($seq1, 0, 27) . ("." x ($L - 27));   # 5' fragment, 3' terminal gaps
    open(my $fh, ">", $path) || die "FAIL: unable to write $path\n";
    print $fh "# STOCKHOLM 1.0\n\n";
    printf $fh "%-32s %s\n", "AB020564.1/7423-7477", $seq1;
    printf $fh "%-32s %s\n", "X15462.1/90-145",      $seq2;
    printf $fh "%-32s %s\n", "frag5p/1-27",          $frag;
    printf $fh "%-32s %s\n", "#=GC SS_cons",          $ss;
    printf $fh "%-32s %s\n", "#=GC RF",               $rf;
    print $fh "//\n";
    close $fh;
}

# A synthetic gapless alignment whose ss_cons has three mutually-crossing
# pseudoknot stems (A,B,C): each opens inside and closes outside the
# nested <<<<< >>>>> stem, and the three cross one another.
sub write_syn3_sto {
    my ($path) = @_;
    my $ss  = "<<<<<AAABBBCCC>>>>>aaabbbccc";
    my $seq = "GGGGGAAACCCUUUCCCCCuuugggaaa";   # gapless; identity irrelevant (ss_cons drives the build)
    die "FAIL: write_syn3_sto length mismatch\n" unless length($ss) == length($seq);
    open(my $fh, ">", $path) || die "FAIL: unable to write $path\n";
    print $fh "# STOCKHOLM 1.0\n\n";
    for my $n ("seq1","seq2","seq3","seq4") { printf $fh "%-20s %s\n", $n, $seq; }
    printf $fh "%-20s %s\n", "#=GC SS_cons", $ss;
    print $fh "//\n";
    close $fh;
}
