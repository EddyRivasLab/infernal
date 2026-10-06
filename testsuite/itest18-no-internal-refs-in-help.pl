#! /usr/bin/perl

# Test that no program's help output (-h, and --devhelp where it
# exists) contains internal development-tracking references, such as
# "brief 123" or "brief 26_0824-085". Those belong in code comments and
# commit messages, not in user-facing output.
#
# Usage:   ./itest18-no-internal-refs-in-help.pl <builddir> <srcdir> <tmpfile prefix>
# Example: ./itest18-no-internal-refs-in-help.pl ..         ..       tmpfoo
#
# EPN, Mon Oct  5 2026 (w/Claude)

BEGIN {
    $builddir  = shift;
    $srcdir    = shift;
    $tmppfx    = shift;
}

@progs    = ("cmalign", "cmbuild", "cmcalibrate", "cmconvert", "cmemit", "cmfetch", "cmpress", "cmscan", "cmsearch", "cmstat");
@devprogs = ("cmbuild", "cmscan", "cmsearch");    # programs with a --devhelp option

foreach $prog (@progs) { if (! -x "$builddir/src/$prog") { die "FAIL: didn't find $prog executable in $builddir/src\n"; } }

$nbad = 0;
foreach $prog (@progs)    { $nbad += check_help($prog, "-h"); }
foreach $prog (@devprogs) { $nbad += check_help($prog, "--devhelp"); }
if ($nbad > 0) { die "FAIL: $nbad help line(s) contain internal tracking references\n"; }

print "ok\n";
exit 0;

# check_help(<prog>, <helpopt>): returns the number of offending lines,
# printing each one.
sub check_help
{
    my ($prog, $helpopt) = @_;
    my @lines;
    my $n = 0;

    @lines = `$builddir/src/$prog $helpopt 2>&1`;
    foreach my $line (@lines) {
        if ($line =~ /\bbriefs?\s*\d/i || $line =~ /\b\d\d_\d{4}\b/) {
            print "$prog $helpopt: $line";
            $n++;
        }
    }
    return $n;
}
