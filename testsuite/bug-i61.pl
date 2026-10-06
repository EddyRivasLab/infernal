#! /usr/bin/perl

# bug i61 - cmalign --mapali corrupts seed rows that begin with '~'
#       (missing data, e.g. fragment ends): residues shifted 5' by
#       the number of leading '~', first residues dropped, and the
#       row truncated so the output is not a valid alignment.
#
# EPN, Tue Oct  6 2026 (w/Claude)
#

$usage = "perl bug-i61.pl <cmbuild> <cmalign>\n";
if ($#ARGV != 1) { die "Wrong argument number.\n$usage"; }

$cmbuild = shift;
$cmalign = shift;
$ok      = 1;

# Make our test alignment, i61.1: row frag1 has 16 leading and 8 trailing '~',
# row frag2 has 2 leading '~'.
#
open (OUT, ">i61.1") || die;
print OUT <<END;
# STOCKHOLM 1.0

frag1     ~~~~~~~~~~~~~~~~UGGG.AGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAG~~~~~~~~
full1     UCCGAUAUAGUGUAAC.GGCUAUCACAUCACGCUUUCACCGUGGAGA.CCGGGGUUCGACUCCCCGUAUCGGAG
full2     UCCGUGAUAGUUUAAU.GGUCAGAAUGGGCGCUUGUCGCGUGCCAGA.UCGGGGUUCAAUUCCCCGUCGCGGAG
frag2     ~~GCACAUGGCGCAGUUGGU.AGCGCGCUUCCCUUGCAAGGAAGAGGUCAUCGGUUCGAUUCCGGUUGCGUCCA
#=GC SS_cons      <<<<<<<..<<<<.........>>>>.<<<<<.......>>>>>.....<<<<<.......>>>>>>>>>>>>.
#=GC RF           [accep]======[=Dloop=]============acd=======[vlp]=====[Tloop]=====[accep]=
//
END
close OUT;
# Make our test sequence to align, i61.2
#
open (OUT, ">i61.2") || die;
print OUT <<END;
>target
GGAGCUAUAGCUCAAUGGCAGAGCGUUUGGCUGACAUCCAAAAGGUUAUGGGUUCGAUUC
END
close OUT;

# build the model
if ($ok) {
  $output = `$cmbuild -F --hand i61.cm i61.1`;
  if ($? != 0) { $ok = 0; }
  # make sure output includes a 'CPU' string indicating success
  if($output !~ /\n\# CPU/) {
    $ok = 0;
  }
}

# align the sequence with --mapali, one-line-per-seq output
if ($ok) {
  $output = `$cmalign --mapali i61.1 --outformat pfam -o i61.3 i61.cm i61.2`;
  if ($? != 0) { $ok = 0; }
  # make sure output includes a 'CPU' string indicating success
  if($output !~ /\n\# CPU/) {
    $ok = 0;
  }
}

# make sure each --mapali row is unchanged apart from '~' and '.' becoming '-'
if ($ok) {
  %expected = ();
  open(IN, "i61.1") || die;
  while($line = <IN>) {
    if($line =~ /^(frag\d|full\d)\s+(\S+)/) {
      ($name, $aseq) = ($1, $2);
      $aseq =~ tr/~./--/;
      $expected{$name} = $aseq;
    }
  }
  close(IN);

  $nfound = 0;
  if(open(IN, "i61.3")) {
    while($line = <IN>) {
      if($line =~ /^(frag\d|full\d)\s+(\S+)/) {
        ($name, $aseq) = ($1, $2);
        $aseq =~ tr/~./--/;
        if($aseq ne $expected{$name}) { $ok = 0; }
        $nfound++;
      }
    }
    close(IN);
  }
  if($nfound != 4) { $ok = 0; }
}

foreach $tmpfile ("i61.1", "i61.2", "i61.3", "i61.cm") {
  unlink $tmpfile if -e $tmpfile;
}
if ($ok) { print "ok\n";     exit 0; }
else     { print "FAILED\n"; exit 1; }
