#! /usr/bin/perl

# bug i57 - cmbuild --eset 0 writes a filter HMM with EFFN 0.000000, which
# HMMER's HMM reader rejects (EFFN must be > 0), so cmsearch/cmscan cannot
# load the resulting .cm file.
#
# This affects any model whose filter HMM is the ML HMM derived from the CM:
# automatically, for zero-basepair models (used here), and also for any
# model built with --p7ml. Basepaired models built without --p7ml are not
# affected, since their filter HMM is built separately and never inherits
# the CM's eff_nseq.
#
# Before the fix: cmbuild succeeds (eff_nseq==0 is not rejected at build
# time), but cmsearch fails to load the model with:
#   Error: bad file format in CM file <f>
#   bad file format for HMM filter
#
# brief 26_0824-079.
#

$usage = "perl bug-i57.pl <cmbuild> <cmsearch>\n";
if ($#ARGV != 1) { die "Wrong argument number.\n$usage"; }

$cmbuild  = shift;
$cmsearch = shift;
$ok       = 1;

# A minimal zero-basepair alignment (no SS_cons pairs).
open (OUT, ">eset0effn.1") || die;
print OUT <<END;
# STOCKHOLM 1.0

seq1               ACGUACGUACGUACGUACGU
seq2               ACGUACGUACGUACGUACGU
#=GC RF            xxxxxxxxxxxxxxxxxxxx
//
END
close OUT;

open (OUT, ">eset0effn.2") || die;
print OUT <<END;
>target1
ACGUACGUACGUACGUACGUACGUACGUACGUACGUACGU
END
close OUT;

if ($ok) {
    system("$cmbuild -F --hand --noss --eset 0 eset0effn.cm eset0effn.1 > /dev/null 2> /dev/null");
    if ($? != 0) { $ok = 0; }
}
if ($ok) {
    # --hmmonly so we don't need to calibrate the model first
    system("$cmsearch --hmmonly eset0effn.cm eset0effn.2 > /dev/null 2> /dev/null");
    if ($? != 0) { $ok = 0; }
}

foreach $tmpfile ("eset0effn.1", "eset0effn.2", "eset0effn.cm") {
    unlink $tmpfile if -e $tmpfile;
}

if ($ok) { print "ok\n";     exit 0; }
else     { print "FAILED\n"; exit 1; }
