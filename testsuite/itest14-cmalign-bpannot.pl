#! /usr/bin/perl

# Test the cmalign --bpstatus / --bpcons structure-status annotation:
#   - one COMBINED #=GR <seq> PS line per sequence (no separate MM line);
#   - a dense #=GC bp_cons family line;
#   - NON-BLANK placeholders so the annotated output is re-readable by the
#     standard Easel Stockholm parser (the write-then-re-read gate, the hole
#     that hid the round-trip bug found in review of the first implementation);
#   - the per-seq PS line survives the low-memory / multi-block (merge) output
#     path, while #=GC bp_cons is suppressed there with a note.
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
                   "bs.sto","bc.sto","def.sto",
                   "tiny.sto","tiny.cm","big.fa","big.sto","big.err","big.sto.rr","log");
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
run("$cmalign --cpu 0 --bpstatus --bpcons $tmppfx.pk.cm $tmppfx.pk.fa > $tmppfx.pk.sto 2>$tmppfx.log",
    "cmalign --bpstatus --bpcons (pknot) failed");

my $nseq_pk = count_match("$tmppfx.pk.fa", qr/^>/);
assert_count("$tmppfx.pk.sto", qr/^#=GR\s+\S+\s+PS\b/, $nseq_pk, "subtest 1: one combined PS line per sequence");
assert_count("$tmppfx.pk.sto", qr/^#=GR\s+\S+\s+MM\b/, 0,        "subtest 1: no separate MM line");
assert_count("$tmppfx.pk.sto", qr/^#=GC\s+bp_cons\b/,  1,        "subtest 1: one #=GC bp_cons line");
assert_dense_ps("$tmppfx.pk.sto", "subtest 1");
# pknot model must show at least one pknot pair-status mark (= / $ / x):
ps_marks_include("$tmppfx.pk.sto", qr/[=\$x]/, "subtest 1: pknot pair-status mark (=/\$/x)");
# bp_cons must carry at least one conservation digit:
gc_has("$tmppfx.pk.sto", "bp_cons", qr/[0-9*]/, "subtest 1: bp_cons conservation digit");
# WRITE-THEN-RE-READ gate (standard Easel reader, via esl-reformat):
reread_ok("$tmppfx.pk.sto", "$tmppfx.pk.rr", "subtest 1 (pknot, both flags)");

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
reread_ok("$tmppfx.bc.sto", "$tmppfx.pk.rr", "subtest 4 (--bpcons only)");

######################################################################
# Subtest 5: default (no flags) output carries NONE of the new lines.
######################################################################
run("$cmalign --cpu 0 $tmppfx.pk.cm $tmppfx.pk.fa > $tmppfx.def.sto 2>$tmppfx.log",
    "cmalign (no flags) failed");
assert_count("$tmppfx.def.sto", qr/^#=GR\s+\S+\s+PS\b/, 0, "subtest 5: no PS lines without flags");
assert_count("$tmppfx.def.sto", qr/^#=GR\s+\S+\s+MM\b/, 0, "subtest 5: no MM lines without flags");
assert_count("$tmppfx.def.sto", qr/^#=GC\s+bp_cons\b/,  0, "subtest 5: no bp_cons without flags");

######################################################################
# Subtest 6: multi-block (merge) path. Emit > 10000 seqs so cmalign routes
# output through the temporary merge file. The per-seq PS line must survive the
# regurgitation and re-read; #=GC bp_cons is suppressed there with a note.
# A tiny model is built just for this volume test (keeps the 10001-seq
# alignment cheap; pseudoknots are not needed to exercise the merge round-trip).
######################################################################
write_tiny_sto("$tmppfx.tiny.sto");
run("$cmbuild --wnone -F $tmppfx.tiny.cm $tmppfx.tiny.sto > $tmppfx.log 2>&1",
    "cmbuild of tiny model failed");
my $nbig = 10001;   # just over CMALIGN_MAX_NSEQ (10000) -> forces a 2nd block / merge path
run("$cmemit -N $nbig --seed 3 -o $tmppfx.big.fa $tmppfx.tiny.cm > $tmppfx.log 2>&1",
    "cmemit of $nbig seqs failed");
run("$cmalign --cpu 0 --bpstatus --bpcons $tmppfx.tiny.cm $tmppfx.big.fa > $tmppfx.big.sto 2>$tmppfx.big.err",
    "cmalign (multi-block, both flags) failed");
assert_count("$tmppfx.big.sto", qr/^#=GR\s+\S+\s+PS\b/, $nbig, "subtest 6: all PS lines survive the merge path");
assert_count("$tmppfx.big.sto", qr/^#=GC\s+bp_cons\b/,  0,     "subtest 6: bp_cons suppressed on the merge path");
slurp("$tmppfx.big.err") =~ /bp_cons.*suppressed|suppressed.*bp_cons/i
    or die "FAIL: subtest 6: expected a bp_cons-suppression note on stderr\n";
reread_ok("$tmppfx.big.sto", "$tmppfx.big.sto.rr", "subtest 6 (multi-block, merge path)");
unlink "$tmppfx.big.sto.rr";

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

# Write-then-re-read gate: the standard Easel Stockholm reader (via esl-reformat)
# must parse the annotated alignment cleanly.
sub reread_ok {
    my ($file, $out, $where) = @_;
    my $err = `$reformat stockholm $file > $out 2>&1`;
    if ($? != 0) { die "FAIL: $where: annotated output is not re-readable by the standard Stockholm parser:\n$err\n"; }
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
