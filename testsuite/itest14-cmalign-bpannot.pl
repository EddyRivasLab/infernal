#! /usr/bin/perl

# Test the cmalign --bpstatus / --bpcons / --bpcov structure-status annotation:
#   - one COMBINED #=GR <seq> PS line per sequence (no separate MM line);
#   - a dense #=GC bp_cons family conservation line;
#   - a dense #=GC bp_cov  family covariation (mutual information) line;
#   - NON-BLANK placeholders so the annotated output is re-readable by the
#     standard Easel Stockholm parser (the write-then-re-read gate, the hole
#     that hid the round-trip bug found in review of the first implementation);
#   - the per-seq PS line and BOTH #=GC family lines survive the low-memory /
#     multi-block (merge) output path (bp_cons / bp_cov via cross-block
#     accumulation, briefs 020 / 022).
#
# Exercises a model that carries BOTH nested and pseudoknot pairs (PK-HAV) and a
# pure-nested model (tRNA), so all of v / ? / : / + / = / $ / x can appear.
#
# Usage:    ./itest14-cmalign-bpannot.pl <builddir> <srcdir> <tmpfile prefix>
# Example:  ./itest14-cmalign-bpannot.pl ..         ..       tmpfoo

use strict;
use warnings;

my $builddir = shift;
my $srcdir   = shift;
my $tmppfx   = shift;

my $cmbuild  = "$builddir/src/cmbuild";
my $cmemit   = "$builddir/src/cmemit";
my $cmalign  = "$builddir/src/cmalign";
my $reformat = "$builddir/easel/miniapps/easel reformat";

if (! -x "$builddir/src/cmbuild")            { die "FAIL: didn't find cmbuild in $builddir/src\n";  }
if (! -x "$builddir/src/cmemit")             { die "FAIL: didn't find cmemit in $builddir/src\n";   }
if (! -x "$builddir/src/cmalign")            { die "FAIL: didn't find cmalign in $builddir/src\n";  }
if (! -x "$builddir/easel/miniapps/easel")   { die "FAIL: didn't find easel in $builddir/easel/miniapps\n"; }
if (! -r "$srcdir/testsuite/PK-HAV.sto")     { die "FAIL: didn't find PK-HAV.sto in $srcdir/testsuite\n"; }
if (! -r "$srcdir/testsuite/tRNA.c.cm")      { die "FAIL: didn't find tRNA.c.cm in $srcdir/testsuite\n";  }

my $pkseed = "$srcdir/testsuite/PK-HAV.sto";
my $trcm   = "$srcdir/testsuite/tRNA.c.cm";

my @tmpsuffixes = ("pk.cm","pk.fa","pk.sto","pk.rr","tr.fa","tr.sto","tr.rr",
                   "bs.sto","bc.sto","bv.sto","def.sto",
                   "tiny.sto","tiny.cm","big.fa","big.sto","big.err","big.sto.rr","log",
                   "sub.fa","sub.sto","sub.rr","miss.fa","miss.sto","miss.rr");
sub cleanup { for my $s (@tmpsuffixes) { unlink "$tmppfx.$s" if -e "$tmppfx.$s"; } }
cleanup();

######################################################################
# Build a pknot CM (PK-HAV carries one A/a pseudoknot stem plus a nested stem).
######################################################################
run("$cmbuild --wnone -F $tmppfx.pk.cm $pkseed > $tmppfx.log 2>&1",
    "cmbuild of PK-HAV pknot seed failed");
file_has_header("$tmppfx.pk.cm", "PKNOT", "yes")
    or die "FAIL: built PK-HAV CM is not a pseudoknot model (PKNOT yes missing)\n";

######################################################################
# Subtest 1: pknot model, single-block, both flags.
######################################################################
run("$cmemit -N 5 --seed 7 -o $tmppfx.pk.fa $tmppfx.pk.cm > $tmppfx.log 2>&1",
    "cmemit from pknot CM failed");
run("$cmalign --cpu 0 --bpstatus --bpcons --bpcov $tmppfx.pk.cm $tmppfx.pk.fa > $tmppfx.pk.sto 2>$tmppfx.log",
    "cmalign --bpstatus --bpcons --bpcov (pknot) failed");

my $nseq_pk = count_match("$tmppfx.pk.fa", qr/^>/);
assert_count("$tmppfx.pk.sto", qr/^#=GR\s+\S+\s+PS\b/, $nseq_pk, "subtest 1: one combined PS line per sequence");
assert_count("$tmppfx.pk.sto", qr/^#=GR\s+\S+\s+MM\b/, 0,        "subtest 1: no separate MM line");
assert_count("$tmppfx.pk.sto", qr/^#=GC\s+bp_cons\b/,  1,        "subtest 1: one #=GC bp_cons line");
assert_count("$tmppfx.pk.sto", qr/^#=GC\s+bp_cov\b/,   1,        "subtest 1: one #=GC bp_cov line");
assert_dense_ps("$tmppfx.pk.sto", "subtest 1");
# both family lines must be dense (whitespace-free token, length == alen):
assert_dense_gc("$tmppfx.pk.sto", "bp_cons", "subtest 1");
assert_dense_gc("$tmppfx.pk.sto", "bp_cov",  "subtest 1");
# pknot model must show at least one pknot pair-status mark (= / $ / x):
ps_marks_include("$tmppfx.pk.sto", qr/[=\$x]/, "subtest 1: pknot pair-status mark (=/\$/x)");
# bp_cons must carry at least one conservation digit; bp_cov a covariation digit:
gc_has("$tmppfx.pk.sto", "bp_cons", qr/[0-9*]/, "subtest 1: bp_cons conservation digit");
gc_has("$tmppfx.pk.sto", "bp_cov",  qr/[0-9*]/, "subtest 1: bp_cov covariation digit");
# bp_cov annotates EXACTLY the same columns as bp_cons (same pair set), but the
# values differ (covariation != conservation) -- a value spot-check:
assert_same_pair_cols("$tmppfx.pk.sto", "bp_cons", "bp_cov", "subtest 1");
gc_lines_differ("$tmppfx.pk.sto", "bp_cons", "bp_cov", "subtest 1: bp_cov differs from bp_cons");
# WRITE-THEN-RE-READ gate (standard Easel reader, via esl-reformat):
reread_ok("$tmppfx.pk.sto", "$tmppfx.pk.rr", "subtest 1 (pknot, all three flags)");

######################################################################
# Subtest 2: pure-nested model (tRNA), single-block, both flags.
######################################################################
run("$cmemit -N 5 --seed 9 -o $tmppfx.tr.fa $trcm > $tmppfx.log 2>&1",
    "cmemit from nested tRNA CM failed");
run("$cmalign --cpu 0 --bpstatus --bpcons $trcm $tmppfx.tr.fa > $tmppfx.tr.sto 2>$tmppfx.log",
    "cmalign --bpstatus --bpcons (nested) failed");
my $nseq_tr = count_match("$tmppfx.tr.fa", qr/^>/);
assert_count("$tmppfx.tr.sto", qr/^#=GR\s+\S+\s+PS\b/, $nseq_tr, "subtest 2: one combined PS line per sequence");
assert_count("$tmppfx.tr.sto", qr/^#=GR\s+\S+\s+MM\b/, 0,        "subtest 2: no separate MM line");
assert_count("$tmppfx.tr.sto", qr/^#=GC\s+bp_cons\b/,  1,        "subtest 2: one #=GC bp_cons line");
assert_dense_ps("$tmppfx.tr.sto", "subtest 2");
reread_ok("$tmppfx.tr.sto", "$tmppfx.tr.rr", "subtest 2 (nested, both flags)");

######################################################################
# Subtest 3: --bpstatus alone (no bp_cons); Subtest 4: --bpcons alone.
######################################################################
run("$cmalign --cpu 0 --bpstatus $tmppfx.pk.cm $tmppfx.pk.fa > $tmppfx.bs.sto 2>$tmppfx.log",
    "cmalign --bpstatus (only) failed");
assert_count("$tmppfx.bs.sto", qr/^#=GR\s+\S+\s+PS\b/, $nseq_pk, "subtest 3: PS lines present with --bpstatus only");
assert_count("$tmppfx.bs.sto", qr/^#=GC\s+bp_cons\b/,  0,        "subtest 3: no bp_cons with --bpstatus only");
reread_ok("$tmppfx.bs.sto", "$tmppfx.pk.rr", "subtest 3 (--bpstatus only)");

run("$cmalign --cpu 0 --bpcons $tmppfx.pk.cm $tmppfx.pk.fa > $tmppfx.bc.sto 2>$tmppfx.log",
    "cmalign --bpcons (only) failed");
assert_count("$tmppfx.bc.sto", qr/^#=GR\s+\S+\s+PS\b/, 0,        "subtest 4: no PS lines with --bpcons only");
assert_count("$tmppfx.bc.sto", qr/^#=GC\s+bp_cons\b/,  1,        "subtest 4: bp_cons present with --bpcons only");
assert_count("$tmppfx.bc.sto", qr/^#=GC\s+bp_cov\b/,   0,        "subtest 4: no bp_cov with --bpcons only");
reread_ok("$tmppfx.bc.sto", "$tmppfx.pk.rr", "subtest 4 (--bpcons only)");

# Subtest 4b: --bpcov alone (independent of --bpcons): bp_cov present, no bp_cons,
# no PS; dense; carries a covariation digit; re-readable.
run("$cmalign --cpu 0 --bpcov $tmppfx.pk.cm $tmppfx.pk.fa > $tmppfx.bv.sto 2>$tmppfx.log",
    "cmalign --bpcov (only) failed");
assert_count("$tmppfx.bv.sto", qr/^#=GR\s+\S+\s+PS\b/, 0,        "subtest 4b: no PS lines with --bpcov only");
assert_count("$tmppfx.bv.sto", qr/^#=GC\s+bp_cons\b/,  0,        "subtest 4b: no bp_cons with --bpcov only");
assert_count("$tmppfx.bv.sto", qr/^#=GC\s+bp_cov\b/,   1,        "subtest 4b: bp_cov present with --bpcov only");
assert_dense_gc("$tmppfx.bv.sto", "bp_cov", "subtest 4b");
gc_has("$tmppfx.bv.sto", "bp_cov", qr/[0-9*]/, "subtest 4b: bp_cov covariation digit");
reread_ok("$tmppfx.bv.sto", "$tmppfx.pk.rr", "subtest 4b (--bpcov only)");

######################################################################
# Subtest 5: default (no flags) output carries NONE of the new lines.
######################################################################
run("$cmalign --cpu 0 $tmppfx.pk.cm $tmppfx.pk.fa > $tmppfx.def.sto 2>$tmppfx.log",
    "cmalign (no flags) failed");
assert_count("$tmppfx.def.sto", qr/^#=GR\s+\S+\s+PS\b/, 0, "subtest 5: no PS lines without flags");
assert_count("$tmppfx.def.sto", qr/^#=GR\s+\S+\s+MM\b/, 0, "subtest 5: no MM lines without flags");
assert_count("$tmppfx.def.sto", qr/^#=GC\s+bp_cons\b/,  0, "subtest 5: no bp_cons without flags");
assert_count("$tmppfx.def.sto", qr/^#=GC\s+bp_cov\b/,   0, "subtest 5: no bp_cov without flags");

######################################################################
# Subtest 6: multi-block (merge) path. Emit > 10000 seqs so cmalign routes
# output through the temporary merge file. The per-seq PS line must survive the
# regurgitation and re-read; #=GC bp_cons is now EMITTED there too (brief 020),
# computed by cross-block accumulation: it must be present, dense (no embedded
# blanks), carry a conservation digit, and re-read cleanly. A tiny model is built
# just for this volume test (keeps the 10001-seq alignment cheap).
######################################################################
write_tiny_sto("$tmppfx.tiny.sto");
run("$cmbuild --wnone -F $tmppfx.tiny.cm $tmppfx.tiny.sto > $tmppfx.log 2>&1",
    "cmbuild of tiny model failed");
my $nbig = 10001;   # just over CMALIGN_MAX_NSEQ (10000) -> forces a 2nd block / merge path
run("$cmemit -N $nbig --seed 3 -o $tmppfx.big.fa $tmppfx.tiny.cm > $tmppfx.log 2>&1",
    "cmemit of $nbig seqs failed");
run("$cmalign --cpu 0 --bpstatus --bpcons --bpcov $tmppfx.tiny.cm $tmppfx.big.fa > $tmppfx.big.sto 2>$tmppfx.big.err",
    "cmalign (multi-block, all three flags) failed");
assert_count("$tmppfx.big.sto", qr/^#=GR\s+\S+\s+PS\b/, $nbig, "subtest 6: all PS lines survive the merge path");
assert_count("$tmppfx.big.sto", qr/^#=GC\s+bp_cons\b/,  1,     "subtest 6: bp_cons emitted on the merge path");
assert_count("$tmppfx.big.sto", qr/^#=GC\s+bp_cov\b/,   1,     "subtest 6: bp_cov emitted on the merge path");
# both family lines must be dense (single whitespace-free token, length == alen):
assert_dense_gc("$tmppfx.big.sto", "bp_cons", "subtest 6");
assert_dense_gc("$tmppfx.big.sto", "bp_cov",  "subtest 6");
# and carry at least one digit (the tiny model's conserved + covarying stems):
gc_has("$tmppfx.big.sto", "bp_cons", qr/[0-9*]/, "subtest 6: multi-block bp_cons conservation digit");
gc_has("$tmppfx.big.sto", "bp_cov",  qr/[0-9*]/, "subtest 6: multi-block bp_cov covariation digit");
# merge-path bp_cov annotates exactly the same columns as bp_cons (same pair set):
assert_same_pair_cols("$tmppfx.big.sto", "bp_cons", "bp_cov", "subtest 6");
# no suppression note should be printed any more:
slurp("$tmppfx.big.err") =~ /bp_cons.*suppressed|suppressed.*bp_cons/i
    and die "FAIL: subtest 6: a stale bp_cons-suppression note was printed on stderr\n";
# WRITE-THEN-RE-READ gate on the merge-path output WITH the bp_cons line present:
reread_ok("$tmppfx.big.sto", "$tmppfx.big.sto.rr", "subtest 6 (multi-block, merge path)");
unlink "$tmppfx.big.sto.rr";

######################################################################
# Subtest 7 (brief 023, Part A): --sub partial-span. Fragmentary targets that
# cover only PART of the model build per-seq sub-CMs that span only a sub-range
# of the consensus. The annotation is computed post-hoc against the FULL cm
# (emap/ct/cmcons), so the marks must stay correct and the line re-readable:
#   - the cpos!=clen RF-length guard must NOT trip on a legitimate partial-span
#     parse (cmalign exits 0);
#   - one dense combined PS line per seq, re-readable by the standard parser;
#   - out-of-span consensus columns carry the no-residue '-' deletion placeholder
#     (NOT a spurious mark), so each fragment's PS line has a long run of '-';
#   - a pair straddling the span boundary (one half in-span, the other out) is a
#     half-present pair -> 'v' (nested) appears.
# 5'-half and 3'-half fragments of the pknot model are the partial-span targets.
######################################################################
write_halves("$tmppfx.pk.fa", "$tmppfx.sub.fa");
run("$cmalign -g --sub --notrunc --cpu 0 --bpstatus --bpcons $tmppfx.pk.cm $tmppfx.sub.fa > $tmppfx.sub.sto 2>$tmppfx.log",
    "subtest 7: cmalign -g --sub --notrunc --bpstatus --bpcons (partial-span) failed (cpos!=clen guard?)");
my $nseq_sub = count_match("$tmppfx.sub.fa", qr/^>/);
assert_count("$tmppfx.sub.sto", qr/^#=GR\s+\S+\s+PS\b/, $nseq_sub, "subtest 7: one combined PS line per partial-span seq");
assert_dense_ps("$tmppfx.sub.sto", "subtest 7");
# out-of-span columns: a long run of the '-' no-residue placeholder must appear:
ps_marks_include("$tmppfx.sub.sto", qr/-{8,}/, "subtest 7: out-of-span '-' placeholder run");
# a boundary-straddling nested pair -> half-present 'v':
ps_marks_include("$tmppfx.sub.sto", qr/v/,     "subtest 7: half-present nested 'v' on a span-boundary pair");
reread_ok("$tmppfx.sub.sto", "$tmppfx.sub.rr", "subtest 7 (--sub partial-span)");

######################################################################
# Subtest 8 (brief 023, Part B): truncated-half '?' fidelity on truncated cmalign
# parses. With --miss, cmalign marks a truncated 5'/3' flank with missing-data
# '~'; a base pair straddling the truncation boundary then has one present half
# and one '~' half. The truncated half must be classified '?' (not 'v', not a
# placeholder, not a spurious '=/$/x'):
#   - NESTED pairs: the present half -> '?' (mode-L/R rule, already correct);
#   - PKNOT pairs: the present half -> '?' too (brief 023 fix in
#     annotate_pknot_pairs_str: a '~' half is "truncated", NOT "broken 'x'").
# So under --miss a '?' must appear AT A PKNOT-LETTER column of SS_cons -- the
# direct regression guard for the fix (pre-fix that column showed 'x').
######################################################################
write_halves("$tmppfx.pk.fa", "$tmppfx.miss.fa");
run("$cmalign --miss --cpu 0 --bpstatus --bpcons $tmppfx.pk.cm $tmppfx.miss.fa > $tmppfx.miss.sto 2>$tmppfx.log",
    "subtest 8: cmalign --miss --bpstatus --bpcons (truncated) failed");
assert_count("$tmppfx.miss.sto", qr/^#=GR\s+\S+\s+PS\b/, count_match("$tmppfx.miss.fa", qr/^>/),
    "subtest 8: one combined PS line per truncated seq");
assert_dense_ps("$tmppfx.miss.sto", "subtest 8");
ps_marks_include("$tmppfx.miss.sto", qr/\?/, "subtest 8: truncated-half '?' mark present");
ps_marks_include("$tmppfx.miss.sto", qr/~/,  "subtest 8: missing-data '~' placeholder present");
# THE FIX: a truncated half of a PKNOT pair must be '?' (not 'x'). Assert some PS
# line carries '?' at a column where SS_cons holds a pknot letter (A-Za-z):
ps_mark_at_pknot_col("$tmppfx.miss.sto", '?', "subtest 8: pknot truncated half -> '?' (brief 023 fix)");
reread_ok("$tmppfx.miss.sto", "$tmppfx.miss.rr", "subtest 8 (--miss truncated)");

######################################################################
print "ok\n";
cleanup();
exit 0;

######################################################################
# Helpers
######################################################################
sub run {
    my ($cmd, $msg) = @_;
    my $out = `$cmd`;
    if ($? != 0) { die "FAIL: $msg\n$out\n"; }
}

sub slurp {
    my ($f) = @_;
    open(my $fh, "<", $f) or return "";
    local $/; my $s = <$fh>; close $fh; return defined $s ? $s : "";
}

sub count_match {
    my ($file, $re) = @_;
    my $n = 0;
    open(my $fh, "<", $file) or die "FAIL: cannot open $file\n";
    while (my $l = <$fh>) { $n++ if $l =~ $re; }
    close $fh;
    return $n;
}

sub assert_count {
    my ($file, $re, $want, $where) = @_;
    my $got = count_match($file, $re);
    if ($got != $want) { die "FAIL: $where: expected $want, got $got (in $file)\n"; }
}

# True if file has a "<tag> <value>" header line whose first value token == value.
sub file_has_header {
    my ($file, $tag, $value) = @_;
    open(my $fh, "<", $file) or die "FAIL: cannot open $file\n";
    while (my $l = <$fh>) {
        chomp $l;
        if ($l =~ /^\Q$tag\E\s+(\S+)/) { close $fh; return ($1 eq $value) ? 1 : 0; }
    }
    close $fh;
    return 0;
}

# Every #=GR PS value must be ONE whitespace-free token (no embedded blanks) and
# exactly as long as the alignment (== the #=GC SS_cons value length). This is
# what makes the annotated output re-readable; a blank in the value would split
# it into fewer-than-4 fields and/or shorten it below alen.
sub assert_dense_ps {
    my ($file, $where) = @_;
    my $alen;
    for my $l (split /\n/, slurp($file)) {
        if ($l =~ /^#=GC\s+SS_cons\s+(\S+)\s*$/) { $alen = length($1); last; }
    }
    defined $alen or die "FAIL: $where: no #=GC SS_cons found to determine alen\n";
    my $checked = 0;
    for my $l (split /\n/, slurp($file)) {
        next unless $l =~ /^#=GR\s+\S+\s+PS\s/;
        my @f = split /\s+/, $l;
        if (scalar(@f) != 4) { die "FAIL: $where: PS line is not a single dense token: '$l'\n"; }
        if (length($f[3]) != $alen) { die "FAIL: $where: PS value length ".length($f[3])." != alen $alen: '$l'\n"; }
        $checked++;
    }
    $checked > 0 or die "FAIL: $where: no PS lines found to dense-check\n";
}

# A #=GC <tag> line must be ONE whitespace-free token exactly as long as the
# alignment (== the #=GC SS_cons value length): dense, no embedded blanks, so it
# round-trips the regurgitator and the standard Stockholm reader.
sub assert_dense_gc {
    my ($file, $tag, $where) = @_;
    my ($alen, $val);
    for my $l (split /\n/, slurp($file)) {
        if ($l =~ /^#=GC\s+SS_cons\s+(\S+)\s*$/)        { $alen = length($1); }
        if ($l =~ /^#=GC\s+\Q$tag\E\s+(\S+)\s*$/)       { $val  = $1; }
    }
    defined $alen or die "FAIL: $where: no #=GC SS_cons found to determine alen\n";
    defined $val  or die "FAIL: $where: no dense (whitespace-free) #=GC $tag line found\n";
    if (length($val) != $alen) { die "FAIL: $where: #=GC $tag length ".length($val)." != alen $alen\n"; }
}

# At least one #=GR PS value contains a char matching $re.
sub ps_marks_include {
    my ($file, $re, $where) = @_;
    for my $l (split /\n/, slurp($file)) {
        next unless $l =~ /^#=GR\s+\S+\s+PS\s+(\S+)\s*$/;
        return 1 if $1 =~ $re;
    }
    die "FAIL: $where: no PS line carried the expected mark\n";
}

# The #=GC <tag> value contains a char matching $re.
sub gc_has {
    my ($file, $tag, $re, $where) = @_;
    for my $l (split /\n/, slurp($file)) {
        next unless $l =~ /^#=GC\s+\Q$tag\E\s+(\S+)\s*$/;
        return 1 if $1 =~ $re;
    }
    die "FAIL: $where: #=GC $tag missing or lacked the expected content\n";
}

# Return the (dense) value string of a #=GC <tag> line, or die.
sub gc_value {
    my ($file, $tag, $where) = @_;
    for my $l (split /\n/, slurp($file)) {
        return $1 if $l =~ /^#=GC\s+\Q$tag\E\s+(\S+)\s*$/;
    }
    die "FAIL: $where: no #=GC $tag line found\n";
}

# Two #=GC lines must mark the SAME set of pair columns: a column is '.' in one
# iff it is '.' in the other. (bp_cov and bp_cons annotate the identical pair set;
# both emit '.' off-pair.) A structural spot-check that does not depend on the
# exact digit values.
sub assert_same_pair_cols {
    my ($file, $tag1, $tag2, $where) = @_;
    my $a = gc_value($file, $tag1, $where);
    my $b = gc_value($file, $tag2, $where);
    length($a) == length($b)
        or die "FAIL: $where: #=GC $tag1 / $tag2 differ in length\n";
    for (my $i = 0; $i < length($a); $i++) {
        my $da = (substr($a,$i,1) eq '.') ? 1 : 0;
        my $db = (substr($b,$i,1) eq '.') ? 1 : 0;
        if ($da != $db) {
            die "FAIL: $where: $tag1 / $tag2 mark different columns at col $i ".
                "('".substr($a,$i,1)."' vs '".substr($b,$i,1)."')\n";
        }
    }
}

# Two #=GC lines must NOT be byte-identical (covariation != conservation).
sub gc_lines_differ {
    my ($file, $tag1, $tag2, $where) = @_;
    my $a = gc_value($file, $tag1, $where);
    my $b = gc_value($file, $tag2, $where);
    $a ne $b or die "FAIL: $where: #=GC $tag1 and $tag2 are byte-identical (expected to differ)\n";
}

# Write-then-re-read gate: the standard Easel Stockholm reader (via esl-reformat)
# must parse the annotated alignment cleanly.
sub reread_ok {
    my ($file, $out, $where) = @_;
    my $err = `$reformat stockholm $file > $out 2>&1`;
    if ($? != 0) { die "FAIL: $where: annotated output is not re-readable by the standard Stockholm parser:\n$err\n"; }
}

# Read a FASTA file and write a new one holding the 5' half and the 3' half of
# each sequence (named <name>_5half / <name>_3half). These fragmentary targets
# cover only part of the model: under --sub they build partial-span sub-CMs;
# under --miss (truncation on) they produce 5'/3'-truncated parses. Determinism
# comes from the fixed-seed cmemit that produced the input.
sub write_halves {
    my ($infa, $outfa) = @_;
    open(my $in, "<", $infa) or die "FAIL: cannot open $infa\n";
    my ($name, %seq, @order);
    while (my $l = <$in>) {
        chomp $l;
        if ($l =~ /^>(\S+)/) { $name = $1; push @order, $name; $seq{$name} = ""; }
        elsif (defined $name) { $l =~ s/\s+//g; $seq{$name} .= $l; }
    }
    close $in;
    open(my $out, ">", $outfa) or die "FAIL: cannot write $outfa\n";
    for my $n (@order) {
        my $s = $seq{$n};
        my $h = int(length($s) / 2);
        next if $h < 4;                       # too short to be a useful half
        print $out ">${n}_5half\n", substr($s, 0, $h), "\n";
        print $out ">${n}_3half\n", substr($s, $h),    "\n";
    }
    close $out;
}

# At least one #=GR PS value carries char <mark> at a column where #=GC SS_cons
# holds a pseudoknot letter (A-Z / a-z). Used to assert the brief-023 fix: a
# truncated half of a PKNOT pair is marked '?' there (pre-fix it was 'x').
sub ps_mark_at_pknot_col {
    my ($file, $mark, $where) = @_;
    my $ss;
    for my $l (split /\n/, slurp($file)) {
        if ($l =~ /^#=GC\s+SS_cons\s+(\S+)\s*$/) { $ss = $1; last; }
    }
    defined $ss or die "FAIL: $where: no #=GC SS_cons found\n";
    my @pkcol = grep { substr($ss,$_,1) =~ /[A-Za-z]/ } (0 .. length($ss)-1);
    @pkcol or die "FAIL: $where: SS_cons carries no pknot letters to check\n";
    for my $l (split /\n/, slurp($file)) {
        next unless $l =~ /^#=GR\s+\S+\s+PS\s+(\S+)\s*$/;
        my $ps = $1;
        for my $c (@pkcol) { return 1 if substr($ps,$c,1) eq $mark; }
    }
    die "FAIL: $where: no PS line carried '$mark' at a pknot column\n";
}

# A tiny gapless nested-stem alignment, used only for the multi-block volume
# test so the 10001-sequence alignment stays cheap.
sub write_tiny_sto {
    my ($path) = @_;
    my $ss  = "<<<<....>>>>..<<<<....>>>>";
    my @sq  = ("GGGGAAAACCCCAAGGGGAAAACCCC",
               "GGGGAAAACCCCAAGGGGAAAACCCC",
               "GGCGAAAACGCCAAGGCGAAAACGCC");
    open(my $fh, ">", $path) or die "FAIL: cannot write $path\n";
    print $fh "# STOCKHOLM 1.0\n\n";
    my $i = 0;
    for my $s (@sq) { $i++; printf $fh "%-14s %s\n", "tiny$i", $s; }
    printf $fh "%-14s %s\n", "#=GC SS_cons", $ss;
    print $fh "//\n";
    close $fh;
}
