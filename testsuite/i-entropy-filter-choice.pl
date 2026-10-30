#! /usr/bin/perl

# The filter HMM and the HMM-only HMM are derived from the CM (brief
# 26_0821-162), so cmbuild's --enone, --ere and --eset (and --rsearch,
# which implies --enone) set their weighting too, and no separate choice
# of the filter HMM's weighting is asked for (as it was before, brief
# 26_0824-082). The options that configured the separately built filter
# HMM are gone; --p7ml is accepted as a no-op.
#
# Checks that the CM weighting options build on their own, that the
# removed options are rejected (writing no CM file), and the HP7
# transform each build records: a2 with the default prior, a1 with
# --rsearch (no default prior to back out), ml with --noss.
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

# [ options, alignment, expect success?, expected HP7 ]
@cases = (
    [ "--enone",                  "efc.1", 1, "a2" ],
    [ "--ere 0.8",                "efc.1", 1, "a2" ],
    [ "--eset 1",                 "efc.1", 1, "a2" ],
    [ "--rsearch $matfile",       "efc.2", 1, "a1" ],
    [ "",                         "efc.1", 1, "a2" ],
    [ "--p7ml",                   "efc.1", 1, "a2" ],
    [ "--eset 1 --p7ml",          "efc.1", 1, "a2" ],
    [ "--noss --eset 0",          "efc.1", 1, "ml" ],
    [ "--p7eent",                 "efc.1", 0, "" ],
    [ "--ere 0.8 --p7ere 0.45",   "efc.1", 0, "" ],
    [ "--p7hemit",                "efc.1", 0, "" ],
    );

foreach $case (@cases) {
    ($opts, $alifile, $expect_ok, $hp7) = @{$case};
    unlink "efc.cm" if -e "efc.cm";
    $output = `$cmbuild -F $opts efc.cm $alifile 2>&1`;
    $rc = $?;
    if ($expect_ok) {
        if ($rc != 0)    { printf("FAIL: cmbuild $opts: expected success, got exit %d\n", $rc >> 8); $ok = 0; }
        if (! -s "efc.cm") { print "FAIL: cmbuild $opts: no CM file written\n"; $ok = 0; }
        elsif (`grep '^HP7' efc.cm` !~ /^HP7\s+$hp7\s*$/) { print "FAIL: cmbuild $opts: HP7 is not $hp7\n"; $ok = 0; }
    }
    else {
        if ($rc == 0)    { print "FAIL: cmbuild $opts: expected failure, got success\n"; $ok = 0; }
        if (-e "efc.cm") { print "FAIL: cmbuild $opts: CM file written despite failure\n"; $ok = 0; }
    }
}

foreach $tmpfile ("efc.1", "efc.2", "efc.cm") {
    unlink $tmpfile if -e $tmpfile;
}

if ($ok) { print "ok\n"; exit 0; }
else     { print "FAILED [entropy filter choice]\n"; exit 1; }
