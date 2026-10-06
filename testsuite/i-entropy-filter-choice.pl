#! /usr/bin/perl

# cmbuild's --enone, --ere and --eset (and --rsearch, which implies
# --enone) set the CM's entropy weighting only; a basepaired model's
# filter HMM keeps its own. cmbuild must refuse these options unless the
# filter HMM's weighting is chosen explicitly (--p7eent, --p7ere or
# --p7ml), or --noss is used (the filter is then the CM's ML HMM).
#
# Checks that the rule rejects/accepts the expected combinations, that a
# rejected run writes no CM file, and that --p7eent is incompatible with
# --p7ml.
#
# brief 26_0824-082.
#

$usage = "perl i-entropy-filter-choice.pl <cmbuild> <RIBOSUM matrix file>\n";
if ($#ARGV != 1) { die "Wrong argument number.\n$usage"; }

$cmbuild = shift;
$matfile = shift;
$ok      = 1;

# A small basepaired alignment.
open (OUT, ">efc.1") || die;
print OUT <<END;
# STOCKHOLM 1.0

seq1               GGGAAACUUCGGUUUCCC
seq2               GGGAAACUACGGUUUCCC
seq3               GGCAAACUUCGGUUUGCC
#=GC SS_cons       <<<<<<......>>>>>>
//
END
close OUT;

# A one-sequence basepaired alignment, for --rsearch.
open (OUT, ">efc.2") || die;
print OUT <<END;
# STOCKHOLM 1.0

seq1               GGGAAACUUCGGUUUCCC
#=GC SS_cons       <<<<<<......>>>>>>
//
END
close OUT;

# [ options, alignment, expect success? ]
@cases = (
    [ "--enone",                  "efc.1", 0 ],
    [ "--ere 0.8",                "efc.1", 0 ],
    [ "--eset 1",                 "efc.1", 0 ],
    [ "--rsearch $matfile",       "efc.2", 0 ],
    [ "--p7eent --p7ml",          "efc.1", 0 ],
    [ "--enone --p7eent --p7ml",  "efc.1", 0 ],
    [ "",                         "efc.1", 1 ],
    [ "--p7eent",                 "efc.1", 1 ],
    [ "--enone --p7eent",         "efc.1", 1 ],
    [ "--ere 0.8 --p7ere 0.45",   "efc.1", 1 ],
    [ "--eset 1 --p7ml",          "efc.1", 1 ],
    [ "--noss --eset 0",          "efc.1", 1 ],
    [ "--rsearch $matfile --p7eent", "efc.2", 1 ],
    );

foreach $case (@cases) {
    ($opts, $alifile, $expect_ok) = @{$case};
    unlink "efc.cm" if -e "efc.cm";
    $output = `$cmbuild -F $opts efc.cm $alifile 2>&1`;
    $rc = $?;
    if ($expect_ok) {
        if ($rc != 0)    { printf("FAIL: cmbuild $opts: expected success, got exit %d\n", $rc >> 8); $ok = 0; }
        if (! -s "efc.cm") { print "FAIL: cmbuild $opts: no CM file written\n"; $ok = 0; }
    }
    else {
        if ($rc == 0)    { print "FAIL: cmbuild $opts: expected failure, got success\n"; $ok = 0; }
        if (-e "efc.cm") { print "FAIL: cmbuild $opts: CM file written despite failure\n"; $ok = 0; }
        # the entropy-flag rule (not the incompatibility check) must name the choices
        if ($opts !~ /--p7ml/ && $output !~ /--p7eent/) { print "FAIL: cmbuild $opts: error message does not mention --p7eent\n"; $ok = 0; }
    }
}

foreach $tmpfile ("efc.1", "efc.2", "efc.cm") {
    unlink $tmpfile if -e $tmpfile;
}

if ($ok) { print "ok\n"; exit 0; }
else     { print "FAILED [entropy filter choice]\n"; exit 1; }
