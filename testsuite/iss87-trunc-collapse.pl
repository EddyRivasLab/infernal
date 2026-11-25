#! /usr/bin/perl

# github issue 87 - cmalign: with tight HMM bands (high tau, e.g. after
#                   tightening to fit --mxsize), a 5'- or 3'-truncated
#                   sequence could be aligned to a single consensus
#                   column (or none) far below its real score.
#
# Each sequence in iss87.fa is a truncated tRNA. Aligned with
# --tau 0.01 --fixedtau, each should span at least 30 consensus
# positions with a positive bit score. Before the fix, all three
# collapsed to a model span of 1 or 0 at about -25 bits.
#
# EPN, Wed Nov 25 2026 (w/Claude)
#

$usage = "perl iss87-trunc-collapse.pl <cmalign> <path to tRNA.c.cm> <path to iss87.fa>\n";
if ($#ARGV != 2) { die "Wrong argument number.\n$usage"; }

$cmalign  = shift;
$cmfile   = shift;
$seqfile  = shift;
$ok       = 1;
$nseq     = 0;

if ($ok) {
  $output = `$cmalign --cpu 0 --tau 0.01 --fixedtau -o iss87.sto $cmfile $seqfile`;
  if ($? != 0) { $ok = 0; }
  # parse the tabular score lines:
  ## idx  seq name          length  cm from    cm to  trunc    bit sc  avg pp  band calc  alignment      total  mem (Mb)
  #    1  tRNA-sample8_t5       59       10       71     5'     20.43   0.879       0.00       0.00       0.00      0.45
  foreach $line (split(/\n/, $output)) {
    if($line !~ m/^\#/ && $line =~ m/\S/) {
      @el_A = split(/\s+/, $line);
      if($el_A[0] eq "") { shift @el_A; }
      if(scalar(@el_A) >= 8) {
        ($cmfrom, $cmto, $bitsc) = ($el_A[3], $el_A[4], $el_A[6]);
        $nseq++;
        if(($cmto - $cmfrom + 1) < 30) { $ok = 0; }
        if($bitsc <= 0.)               { $ok = 0; }
      }
    }
  }
  if($nseq != 3) { $ok = 0; }
}

foreach $tmpfile ("iss87.sto") {
  unlink $tmpfile if -e $tmpfile;
}

if ($ok) { print "ok\n";     exit 0; }
else     { print "FAILED\n"; exit 1; }
