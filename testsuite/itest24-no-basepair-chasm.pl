#! /usr/bin/perl

# itest24: a CM with no basepairs is built the same way as a CM with
# basepairs, minus the pairs, so adding or removing basepairs does not
# change the rest of the model.
#
# Before brief 26_0821-182, cmbuild built a 0-bp CM as if it were an
# HMM (an HMMER3-like prior, parse-tree doctoring, a 0.38 bit entropy
# weighting target, no separate HMM-only HMM), so adding one basepair
# re-parameterized every column.
#
# From tRNA.sto we build three CMs: S, with its structure; B1, with
# only its first basepair; B0, with none (--noss). Then:
#
# 1. With entropy weighting off (--enone, so the effective sequence
#    number is the same for all three), the parameters of the columns
#    that are unpaired in all three are identical, exactly (to the
#    precision of the CM file):
#    - the match emission scores of each such column (an ML state of a
#      MATL node, or an MR state of a MATR node), and
#    - the transition scores (ML, D and IL states) of each MATL node
#      whose next node is a MATL node for the next column, in all
#      three CMs. Other transitions can't match: their prior depends
#      on the type of the next node.
# 2. With default entropy weighting, each CM gets its own effective
#    sequence number (EFFN), because the paired columns add relative
#    entropy. That is the only difference: B0 built with --eset set to
#    S's EFFN has the match emission scores of S's default build in
#    every column unpaired in S, to within 0.001 (EFFN is written with
#    6 decimals).
# 3. B0 has an HMM-only HMM (HP7 h3pri), as S does, and --noh3pri,
#    which used to turn the HMM-like build off, changes nothing.
#
# Usage:    ./itest24-no-basepair-chasm.pl <builddir> <srcdir> <tmpfile prefix>
# Example:  ./itest24-no-basepair-chasm.pl ..         ..       tmpfoo
#
# EPN, Thu Nov  5 2026 (w/Claude)
# brief 26_0821-182

BEGIN {
    $builddir = shift;
    $srcdir   = shift;
    $tmppfx   = shift;
}

use strict;
use warnings;

our ($builddir, $srcdir, $tmppfx);

my $cmbuild = "$builddir/src/cmbuild";
my $ali     = "$srcdir/testsuite/tRNA.sto";
if (! -x $cmbuild) { die "FAIL: didn't find cmbuild binary in $builddir/src\n"; }
if (! -e $ali)     { die "FAIL: didn't find $ali\n"; }

my @tmpfiles = ();
sub tmp { my $f = "$tmppfx.$_[0]"; push @tmpfiles, $f; unlink $f if -e $f; return $f; }

# B1's alignment: tRNA.sto with all basepairs but the first removed from SS_cons
my $ali1 = tmp("b1.sto");
write_lines($ali1, one_pair(read_lines($ali)));

sub build {
    my ($tag, $opts, $alifile) = @_;
    my $cm = tmp("$tag.cm");
    `$cmbuild --cpu 1 -F $opts $cm $alifile > /dev/null 2>&1`;
    if ($? != 0) { die "FAIL: cmbuild $opts ($tag) failed\n"; }
    return $cm;
}

# 1.
my %cm;
$cm{S}  = build("S.enone",  "--enone",        $ali);
$cm{B1} = build("B1.enone", "--enone",        $ali1);
$cm{B0} = build("B0.enone", "--enone --noss", $ali);
my %p = map { $_ => parse_cm($cm{$_}) } keys %cm;
if (npairs($p{S}) < 2 || npairs($p{B1}) != 1 || npairs($p{B0}) != 0) { die "FAIL: subtest 1: unexpected numbers of basepairs\n"; }
my ($nemit, $ntrans) = (0, 0);
foreach my $col (sort { $a <=> $b } keys %{$p{B0}{emit}}) {
    next if (! exists $p{S}{emit}{$col} || ! exists $p{B1}{emit}{$col});
    foreach my $m ("S", "B1") {
        if ($p{$m}{emit}{$col} ne $p{B0}{emit}{$col}) { die "FAIL: subtest 1: column $col emission scores differ, $m vs B0:\n$p{$m}{emit}{$col}\n$p{B0}{emit}{$col}\n"; }
    }
    $nemit++;
}
foreach my $col (sort { $a <=> $b } keys %{$p{B0}{trans}}) {
    next if (! exists $p{S}{trans}{$col} || ! exists $p{B1}{trans}{$col});
    foreach my $m ("S", "B1") {
        if ($p{$m}{trans}{$col} ne $p{B0}{trans}{$col}) { die "FAIL: subtest 1: column $col transition scores differ, $m vs B0:\n$p{$m}{trans}{$col}\n$p{B0}{trans}{$col}\n"; }
    }
    $ntrans++;
}
# every column unpaired in S is compared (tRNA.sto: 29 of 71)
if ($nemit != scalar(keys %{$p{S}{emit}}) || $nemit < 25 || $ntrans < 10) { die "FAIL: subtest 1: too few columns compared ($nemit emission, $ntrans transition)\n"; }

# 2.
my $sdef  = build("S.def", "", $ali);
my $effn  = header_value($sdef, "EFFN");
if (! defined $effn) { die "FAIL: subtest 2: no EFFN line\n"; }
my $b0set = build("B0.eset", "--noss --eset $effn", $ali);
my $b0def = build("B0.def",  "--noss", $ali);
if (header_value($b0def, "EFFN") == $effn) { die "FAIL: subtest 2: S and B0 have the same default EFFN; the test needs them to differ\n"; }
my $ps    = parse_cm($sdef);
my $pb    = parse_cm($b0set);
my $nset  = 0;
foreach my $col (keys %{$ps->{emit}}) {
    if (! exists $pb->{emit}{$col}) { die "FAIL: subtest 2: column $col unpaired in S but not in B0\n"; }
    my @a = split(" ", $ps->{emit}{$col});
    my @b = split(" ", $pb->{emit}{$col});
    for (my $i = 0; $i < scalar(@a); $i++) {
        if (abs($a[$i] - $b[$i]) > 0.0011) { die "FAIL: subtest 2: column $col emission scores differ beyond 0.001: $ps->{emit}{$col} vs $pb->{emit}{$col}\n"; }
    }
    $nset++;
}
if ($nset < 25) { die "FAIL: subtest 2: too few columns compared ($nset)\n"; }

# 3.
if ((header_value($b0def, "HP7") // "") ne "h3pri") { die "FAIL: subtest 3: B0's HP7 is not h3pri\n"; }
if ((header_value($sdef,  "HP7") // "") ne "h3pri") { die "FAIL: subtest 3: S's HP7 is not h3pri\n"; }
my $b0noh = build("B0.noh3pri", "--noss --noh3pri", $ali);
if (stripped($b0noh) ne stripped($b0def)) { die "FAIL: subtest 3: --noh3pri changed the 0-bp build\n"; }

foreach my $f (@tmpfiles) { unlink $f if -e $f; }
print "ok\n";
exit 0;

sub read_lines  { my ($f) = @_; open(my $fh, "<", $f) || die "FAIL: couldn't open $f\n"; my @l = <$fh>; close $fh; return @l; }
sub write_lines { my ($f, @l) = @_; open(my $fh, ">", $f) || die "FAIL: couldn't open $f\n"; print $fh @l; close $fh; }
sub stripped    { return join("", grep { ! /^(DATE|COM) / } read_lines($_[0])); }
sub header_value {
    my ($f, $tag) = @_;
    foreach my $l (read_lines($f)) { if ($l =~ /^\Q$tag\E\s+(.*?)\s*$/) { return $1; } if ($l =~ /^CM\s*$/) { last; } }
    return undef;
}

# one_pair(<Stockholm lines>): keep only the first basepair in SS_cons
# (its first opening bracket and the bracket that closes it); every
# other bracket becomes '.'
sub one_pair {
    my @lines = @_;
    my $ss = join("", map { /^#=GC SS_cons\s+(\S+)/ ? $1 : () } @lines);
    my @c  = split(//, $ss);
    my ($open, $close, $depth) = (-1, -1, 0);
    for (my $i = 0; $i < scalar(@c); $i++) {
        if ($c[$i] =~ /[<(\[{]/) { if ($open < 0) { $open = $i; } $depth++; }
        elsif ($c[$i] =~ /[>)\]}]/) { $depth--; if ($open >= 0 && $depth == 0) { $close = $i; last; } }
    }
    if ($open < 0 || $close < 0) { die "FAIL: no basepair in $ali\n"; }
    my $pos = 0;
    foreach my $l (@lines) {
        if ($l =~ /^(#=GC SS_cons\s+)(\S+)(.*)$/s) {
            my ($pfx, $s, $sfx) = ($1, $2, $3);
            my @s = split(//, $s);
            for (my $i = 0; $i < scalar(@s); $i++, $pos++) {
                if ($s[$i] =~ /[<>()\[\]{}]/ && $pos != $open && $pos != $close) { $s[$i] = "."; }
            }
            $l = $pfx . join("", @s) . $sfx;
        }
    }
    return @lines;
}

# parse_cm(<cmfile>): for the first CM in the file, return a hash ref:
#   {npairs}:     number of MATP nodes
#   {emit}{col}:  match emission scores of unpaired consensus column col
#                 (ML of a MATL node, MR of a MATR node), as a string
#   {trans}{col}: transition scores of the ML, D and IL states of the
#                 MATL node for column col, if its next node is a MATL
#                 node for column col+1
# A state line is: type, v, plast, pnum, cfirst, cnum, four band
# values (which depend on the whole model), then transition scores
# (cnum of them), then emission scores.
sub parse_cm {
    my ($f) = @_;
    my %r = (npairs => 0, emit => {}, trans => {});
    my @nodes = ();  # [ type, lpos, rpos, { state type => [ trans, emit ] } ]
    my $in = 0;
    foreach my $l (read_lines($f)) {
        if ($l =~ /^CM\s*$/) { $in = 1; next; }
        last if ($in && $l =~ /^\/\//);
        next if (! $in);
        if ($l =~ /^\s+\[\s+(\S+)\s+\d+\s+\]\s+(\S+)\s+(\S+)/) { push @nodes, [ $1, $2, $3, {} ]; next; }
        my @t = split(" ", $l);
        next if (scalar(@t) < 10 || ! @nodes || $t[0] eq "B"); # a B state's cfirst/cnum are its children
        my $cnum  = $t[5];
        my @sc    = @t[10..$#t];
        my @trans = @sc[0..($cnum-1)];
        my @emit  = @sc[$cnum..$#sc];
        $nodes[-1][3]{$t[0]} = [ join(" ", @trans), join(" ", @emit) ];
    }
    for (my $n = 0; $n < scalar(@nodes); $n++) {
        my ($type, $lpos, $rpos, $st) = @{$nodes[$n]};
        if    ($type eq "MATP") { $r{npairs}++; }
        elsif ($type eq "MATL") { $r{emit}{$lpos} = $st->{ML}[1]; }
        elsif ($type eq "MATR") { $r{emit}{$rpos} = $st->{MR}[1]; }
        if ($type eq "MATL" && $n+1 < scalar(@nodes) && $nodes[$n+1][0] eq "MATL" && $nodes[$n+1][1] == $lpos+1) {
            $r{trans}{$lpos} = join(" | ", map { $st->{$_}[0] } ("ML", "D", "IL"));
        }
    }
    return \%r;
}
sub npairs { return $_[0]{npairs}; }
