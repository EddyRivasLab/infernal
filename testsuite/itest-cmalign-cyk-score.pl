#! /usr/bin/perl

# Check that the score cmalign reports for a CYK alignment is the score
# of the parse it reports, for parses that use local ends (EL states)
# and truncated (marginal) alignment modes.
#
# cmalign centers the CM DP scores: it lowers every residue's emission
# score by a constant c for the DP and adds c*L back to the score.
# That is exact only if every residue of every parse gets the shift,
# including residues emitted by EL states and by MP states in marginal
# (Left/Right) modes. The --tfile output scores each parse with
# ParsetreeScore() from the CM's own (unshifted) scores, so a residue
# that missed the shift shows up as a difference of c bits or more
# between the reported score and the --tfile score.
#
# The sequences are three emitted tRNAs (emitted-tRNA.fa), each also
# truncated at the 5' end, the 3' end and both ends, and with a 40 nt
# insertion. The test fails if the local-mode parses don't use EL
# states, or don't use marginal modes, so it can't pass by no longer
# exercising them.
#
# Usage:    ./itest-cmalign-cyk-score.pl <builddir> <srcdir> <tmpfile prefix>
# Example:  ./itest-cmalign-cyk-score.pl ..         ..       tmpfoo
#
# EPN, Fri Oct 30 2026 (w/Claude)

BEGIN {
    $builddir = shift;
    $srcdir   = shift;
    $tmppfx   = shift;
}
use lib "$srcdir/testsuite";  # The BEGIN is necessary to make this work: sets $srcdir at compile-time
use i1;

$tol = 0.015; # both scores are printed with 2 decimals

if (! -x "$builddir/src/cmalign") { die "FAIL: didn't find cmalign binary in $builddir/src\n"; }
$cmfile  = "$srcdir/testsuite/tRNA1415G.c.cm";
$seqfile = "$srcdir/testsuite/emitted-tRNA.fa";
if (! -e $cmfile)  { die "FAIL: didn't find $cmfile\n"; }
if (! -e $seqfile) { die "FAIL: didn't find $seqfile\n"; }

# Create the test sequences
$ins = "GAUCCAAGUUCAGGACUCAAGGUCCAUUGACCGAAUCGAA";
open(SQ, $seqfile) || die "FAIL: couldn't open $seqfile\n";
$name = ""; %seqH = (); @nameA = ();
while ($line = <SQ>) {
    chomp $line;
    if ($line =~ /^>(\S+)/) { $name = $1; push(@nameA, $name); $seqH{$name} = ""; }
    else                    { $seqH{$name} .= $line; }
}
close(SQ);
open(OUT, ">$tmppfx.fa") || die "FAIL: couldn't open $tmppfx.fa for write\n";
foreach $name (@nameA) {
    $s = $seqH{$name};
    $L = length($s);
    printf OUT (">%s\n%s\n",      $name, $s);
    printf OUT (">%s-t5\n%s\n",   $name, substr($s, 25));
    printf OUT (">%s-t3\n%s\n",   $name, substr($s, 0, $L-25));
    printf OUT (">%s-t53\n%s\n",  $name, substr($s, 15, $L-30));
    printf OUT (">%s-ins\n%s\n",  $name, substr($s, 0, 35) . $ins . substr($s, 35));
}
close(OUT);

# Align in local (truncated) mode, local --notrunc mode and global mode
$nEL_local = $nmarg_local = $nEL_notrunc = 0;
foreach $opts ("", "--notrunc", "-g") {
    $output = `$builddir/src/cmalign --cyk --cpu 0 $opts --tfile $tmppfx.tfile -o $tmppfx.stk $cmfile $tmppfx.fa 2>&1`;
    if ($? != 0)               { die "FAIL: cmalign --cyk $opts failed\n"; }
    if ($output !~ /\n\# CPU/) { die "FAIL: cmalign --cyk $opts output incomplete\n"; }

    # reported scores
    %repH = ();
    foreach $line (split(/\n/, $output)) {
        if ($line =~ /^\#/ || $line !~ /\S/) { next; }
        @fA = split(/\s+/, $line);
        if ($fA[0] eq "") { shift @fA; }
        if (scalar(@fA) > 7) { $repH{$fA[1]} = $fA[6]; }
    }

    # parse scores, EL residues and marginal modes from the --tfile
    %tfH = (); %elH = (); %margH = (); $name = "";
    open(TF, "$tmppfx.tfile") || die "FAIL: couldn't open $tmppfx.tfile\n";
    while ($line = <TF>) {
        if    ($line =~ /^>(\S+)/)               { $name = $1; $elH{$name} = 0; $margH{$name} = 0; }
        elsif ($line =~ /^\s+SCORE:\s+(\S+) bits/) { $tfH{$name} = $1; }
        elsif ($line =~ /^\s+\d+\s+(\d+)\S*\s+(\d+)\S*\s+\d+(\S+)\s+(\S+)/) {
            ($i, $j, $st, $mode) = ($1, $2, $3, $4);
            if ($st eq "EL" && $j >= $i)              { $elH{$name} += $j - $i + 1; }
            if ($mode eq "Left" || $mode eq "Right")  { $margH{$name}++; }
        }
    }
    close(TF);

    foreach $name (keys %repH) {
        if (! exists $tfH{$name}) { die "FAIL: no --tfile score for $name (cmalign --cyk $opts)\n"; }
        if (abs($repH{$name} - $tfH{$name}) > $tol) {
            die sprintf("FAIL: cmalign --cyk %s: %s reported score %s, parse score %s\n", $opts, $name, $repH{$name}, $tfH{$name});
        }
        if ($opts eq "")          { if ($elH{$name} > 0) { $nEL_local++;   } if ($margH{$name} > 0) { $nmarg_local++; } }
        if ($opts eq "--notrunc") { if ($elH{$name} > 0) { $nEL_notrunc++; } }
    }
    if (scalar(keys %repH) != 5 * scalar(@nameA)) { die "FAIL: cmalign --cyk $opts aligned the wrong number of sequences\n"; }
}

# the test must exercise EL residues and marginal modes
if ($nEL_local   == 0) { die "FAIL: no local-mode parse uses an EL state\n"; }
if ($nmarg_local == 0) { die "FAIL: no local-mode parse uses a marginal mode\n"; }
if ($nEL_notrunc == 0) { die "FAIL: no --notrunc parse uses an EL state\n"; }

foreach $tmpfile ("$tmppfx.fa", "$tmppfx.tfile", "$tmppfx.stk") { unlink $tmpfile if -e $tmpfile; }

print "ok\n";
exit 0;
