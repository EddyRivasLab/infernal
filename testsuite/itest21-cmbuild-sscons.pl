#! /usr/bin/perl

# Test cmbuild --sscons <tag>: build the consensus structure from the
# #=GC <tag> line instead of #=GC SS_cons.
#
#   - a CM built with --sscons ALT from an alignment whose SS_cons is
#     structure X and whose #=GC ALT is structure Y is the same CM as one
#     built from a copy of the alignment whose SS_cons IS structure Y;
#   - it differs from the default (SS_cons) build;
#   - '--sscons SS_cons' is a no-op;
#   - a missing tag fails with an error naming the tag;
#   - a tag of the wrong length fails (the Stockholm parser itself rejects
#     it, before cmbuild's own length check);
#   - --sscons is incompatible with --noss.
#
# Usage:    ./itest21-cmbuild-sscons.pl <builddir> <srcdir> <tmpfile prefix>
# Example:  ./itest21-cmbuild-sscons.pl ..         ..       tmpfoo

use strict;
use warnings;

my $builddir = shift;
my $srcdir   = shift;
my $tmppfx   = shift;

my $cmbuild  = "$builddir/src/cmbuild";
if (! -x $cmbuild) { die "FAIL: didn't find cmbuild in $builddir/src\n"; }

my @tmpsuffixes = ("alt.sto", "ref.sto", "short.sto", "alt.cm", "ref.cm", "def.cm", "noop.cm", "x.cm", "log");
sub cleanup { for my $s (@tmpsuffixes) { unlink "$tmppfx.$s" if -e "$tmppfx.$s"; } }
cleanup();

#            two hairpins                      one long hairpin
my $ssX  = "<<<<....>>>>..<<<<....>>>>";
my $ssY  = "<<<<<<.............>>>>>>.";
my @sq   = ("GGGGAAAACCCCAAGGGGAAAACCCC",
            "GGGGAAAACCCCAAGGGGAAAACCCC",
            "GGCGAAAACGCCAAGGCGAAAACGCC");

write_sto("$tmppfx.alt.sto",   $ssX, [ ["ALT", $ssY] ]);
write_sto("$tmppfx.ref.sto",   $ssY, [ ]);
write_sto("$tmppfx.short.sto", $ssX, [ ["ALT", substr($ssY, 1)] ]);

# --sscons ALT == building from SS_cons = Y
run("$cmbuild --wnone --sscons ALT -F $tmppfx.alt.cm $tmppfx.alt.sto > $tmppfx.log 2>&1", "cmbuild --sscons ALT failed");
run("$cmbuild --wnone -F $tmppfx.ref.cm $tmppfx.ref.sto > $tmppfx.log 2>&1", "cmbuild of reference alignment failed");
run("$cmbuild --wnone -F $tmppfx.def.cm $tmppfx.alt.sto > $tmppfx.log 2>&1", "cmbuild (default) failed");
run("$cmbuild --wnone --sscons SS_cons -F $tmppfx.noop.cm $tmppfx.alt.sto > $tmppfx.log 2>&1", "cmbuild --sscons SS_cons failed");

my $alt  = model_body("$tmppfx.alt.cm");
my $ref  = model_body("$tmppfx.ref.cm");
my $def  = model_body("$tmppfx.def.cm");
my $noop = model_body("$tmppfx.noop.cm");
if ($alt ne $ref)  { die "FAIL: --sscons ALT model differs from model built with SS_cons = ALT structure\n"; }
if ($alt eq $def)  { die "FAIL: --sscons ALT model is identical to the default (SS_cons) model\n"; }
if ($noop ne $def) { die "FAIL: --sscons SS_cons is not a no-op\n"; }

# errors
my $out;
$out = `$cmbuild --wnone --sscons NOSUCHTAG -F $tmppfx.x.cm $tmppfx.alt.sto 2>&1`;
if ($? == 0)                     { die "FAIL: --sscons with a missing tag did not fail\n"; }
if ($out !~ /NOSUCHTAG/)         { die "FAIL: missing-tag error does not name the tag:\n$out\n"; }
$out = `$cmbuild --wnone --sscons ALT -F $tmppfx.x.cm $tmppfx.short.sto 2>&1`;
if ($? == 0)                     { die "FAIL: --sscons with a wrong-length tag did not fail\n"; }
$out = `$cmbuild --wnone --sscons ALT --noss -F $tmppfx.x.cm $tmppfx.alt.sto 2>&1`;
if ($? == 0)                     { die "FAIL: --sscons with --noss did not fail\n"; }

print "ok\n";
cleanup();
exit 0;

sub run {
    my ($cmd, $msg) = @_;
    my $o = `$cmd`;
    if ($? != 0) { die "FAIL: $msg\n$o\n"; }
}

# Model file contents minus lines that legitimately differ between builds
# (version/date, command line, input file, timing, checksum).
sub model_body {
    my ($f) = @_;
    open(my $fh, "<", $f) or die "FAIL: cannot open $f\n";
    my $s = "";
    while (my $l = <$fh>) {
        next if $l =~ /^(INFERNAL|HMMER|COM|DATE|CKSUM|NAME|EFFN|NSEQ)\b/;
        $s .= $l;
    }
    close $fh;
    return $s;
}

sub write_sto {
    my ($path, $ss, $gcref) = @_;
    open(my $fh, ">", $path) or die "FAIL: cannot write $path\n";
    print $fh "# STOCKHOLM 1.0\n\n";
    my $i = 0;
    for my $s (@sq) { $i++; printf $fh "%-14s %s\n", "seq$i", $s; }
    printf $fh "%-14s %s\n", "#=GC SS_cons", $ss;
    for my $gc (@$gcref) { printf $fh "%-14s %s\n", "#=GC $gc->[0]", $gc->[1]; }
    print $fh "//\n";
    close $fh;
}
