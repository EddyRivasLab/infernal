#! /usr/bin/perl

# Test the cmalign --mapali option.
# 
# Usage:    ./itest1-cmalign-mapali.pl  <builddir> <srcdir> <tmpfile prefix>
# Example:  ./itest1-cmalign-mapali.pl  ..         ..       foo
#
# EPN, Thu Apr 26 14:47:10 2012
# Based on HMMER3's i6-hmmalign-mapali [SRE, Mon May 25 09:52:48 2009]

$builddir  = shift;
$srcdir    = shift;
$tmppfx    = shift;

if (! -x "$builddir/src/cmalign")           { die "FAIL: didn't find cmalign in $builddir/src/";          }
if (! -x "$builddir/easel/miniapps/easel")  { die "FAIL: didn't find easel in $builddir/easel/miniapps/"; }
if (! -r "$srcdir/testsuite/Vault.c.cm")    { die "FAIL: didn't find Vault.c.cm in $srcdir/testsuite/";   }
if (! -r "$srcdir/testsuite/Vault.sto")     { die "FAIL: didn't find Vault.sto in $srcdir/testsuite/";    }

$cmalign      = "$builddir/src/cmalign";
$eslreformat  = "$builddir/easel/miniapps/easel reformat";
$testsuitedir = "$srcdir/testsuite";

system("$eslreformat -u --rename foo fasta $testsuitedir/Vault.sto > $tmppfx.fa");
if ($? != 0)   { print "FAIL: esl-reformat failed unexpectedly\n"; exit 1; }

system("$cmalign -o $tmppfx.sto --mapali $testsuitedir/Vault.sto $testsuitedir/Vault.c.cm $tmppfx.fa");
if ($? != 0)   { print "FAIL: cmalign failed unexpectedly\n"; exit 1; }

system("$eslreformat -u fasta $tmppfx.sto > $tmppfx.2.fa");
if ($? != 0)   { print "FAIL: esl-reformat failed unexpectedly\n"; exit 1; }

system("$eslreformat -u fasta $testsuitedir/Vault.sto > $tmppfx.3.fa");
if ($? != 0)   { print "FAIL: esl-reformat failed unexpectedly\n"; exit 1; }

system("cat $tmppfx.fa >> $tmppfx.3.fa");
if ($? != 0)   { print "FAIL: cat failed unexpectedly\n"; exit 1; }

system("diff $tmppfx.2.fa $tmppfx.3.fa");
if ($? != 0)   { print "FAIL: --mapali doesn't produce expected sequences\n"; exit 1; }

print "ok\n"; 
unlink "$tmppfx.fa";
unlink "$tmppfx.sto";
unlink "$tmppfx.2.fa";
unlink "$tmppfx.3.fa";
exit 0;
