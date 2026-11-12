#! /usr/bin/perl

# bug i47 - cmbuild --p7ml option does nothing
#
# Since the filter HMM became the CM's ML HMM for every model (brief
# 26_0821-162), --p7ml is accepted as a no-op: this test now checks that
# the --p7ml and default models are identical apart from COM and DATE.
#
# EPN, Sat Mar  4 09:38:03 2017
#

$usage = "perl bug-i47.pl <cmbuild>\n";
if ($#ARGV != 0) { die "Wrong argument number.\n$usage"; }

$cmbuild = shift;
$ok      = 1;

# Make our test alignment, i44.1
#
open (OUT, ">i47.1") || die;
print OUT <<END;
# STOCKHOLM 1.0

seq1         AAAACCCCGGGGUUUU
seq2         AAAACCCCGGGGUUUU
seq3         AAAACCCCGGGGUUUU
seq4         AAAACCCCGGGGUUUU
seq5         AAAACCCCGGGGUUUU
seq6         AAAACCCCGGGGUUUU
seq7         AAAACCCCGGGGUUUU
seq8         AAAACCCCGGGGUUUU
seq9         AAAACCCCGGGGUUUU
seq10        AAAACCCCGGGGUUUU
#=GC SS_cons <<..<<<..>>>..>>
#=GC RF      AAAACCCCGGGGUUUU
//
END
close OUT;

# build the default model
if ($ok) { 
  $output = `$cmbuild -F i47.1.cm i47.1`;
  if ($? != 0) { $ok = 0; }
  # make sure output includes a 'CPU' string indicating success
  if($output !~ /\n\# CPU/) { 
    $ok = 0;
  }
  # remove COM lines
  if($ok) { 
    `grep -v -e ^COM -e ^DATE i47.1.cm > i47.3.cm`;
    if ($? != 0) { $ok = 0; }
  }
}
# build the --p7ml model
if ($ok) { 
  $output = `$cmbuild -F --p7ml i47.2.cm i47.1`;
  if ($? != 0) { $ok = 0; }
  # make sure output includes a 'CPU' string indicating success
  if($output !~ /\n\# CPU/) { 
    $ok = 0;
  }
  # remove COM lines
  if($ok) { 
    `grep -v -e ^COM -e ^DATE i47.2.cm > i47.4.cm`; 
    if ($? != 0) { $ok = 0; }
  }
}

if($ok) { 
  # --p7ml is a no-op: both files must be equal and diff will return nothing
  $ndifflines = `diff i47.3.cm i47.4.cm | wc -l`;
  if ($? != 0) { $ok = 0; }
  if($ok) { 
    $ok = ($ndifflines == 0) ? 1 : 0;
  }
}

foreach $tmpfile ("i47.1", "i47.1.cm", "i47.2.cm", "i47.3.cm", "i47.4.cm") { 
  unlink $tmpfile if -e $tmpfile;
}
if ($ok) { print "ok\n";     exit 0; }
else     { print "FAILED\n"; exit 1; }


