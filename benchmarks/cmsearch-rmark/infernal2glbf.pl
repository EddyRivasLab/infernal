#! /usr/bin/perl

# @LICENSE@
# EPN 09.06.05
# modified from hmmer2gdf from SRE

# Usage: perl infernal2glbf.pl <cmsearch output file>
#
# Converts cmsearch output to GLBF.
# GLBF format is just like GLF format but with bounds of hits
# <seq name> <score> <start posn> <end posn> <orientation 0 (forward) of 1 (reverse)>
# Score is either E-value (if cmsearch was run with stats enabled), else it's bit score 
# Order is sorted by sequence position. 
#
#
# Options from :
#    -E <x>         : use E values [default], sets max E-val to keep as <x> [default = 2]
#    -B <x>         : use bit scores, sets min score to keep as <x>
#
# SRE, Wed Oct 28 14:03:52 1998
# RCS $Id$
#

use Getopt::Std;
use infernal;

$use_evalues   = 1;
$use_bitscores = 0;
$e_cutoff =   2;
$b_cutoff = 0.0;

getopts('E:B:');
if (defined $opt_E) { $e_cutoff = $opt_E; }
if (defined $opt_B) { $b_cutoff = $opt_B; $use_evalues = 0; $use_bitscores = 1; }

$output = join("",<>);
&infernal::ParseINFERNAL($output);

# First determine if infernal was run with or without E-values
if($infernal::nhit > 0)
{
    if (exists($infernal::hitevalue[0]))
    {
	$has_evalues = 1;
    }
    else
    {
	$has_evalues = 0;
    }
}

if($use_evalues && (!$has_evalues))
{
    die("ERROR in infernal2glbf.pl, trying to use E-values but none reported.\n");
}

for ($i = 0; $i < $infernal::nhit; $i++)
{
    if ((($use_bitscores) && $infernal::hitbitscore[$i] > $b_cutoff) ||
	(($use_evalues)  && $infernal::hitevalue[$i]   < $e_cutoff))
    {
	#printf("%-24s %-6f\n", $infernal::targname[$i], $infernal::seqbitscore{$infernal::targname[$i]}); 
	if($infernal::hitsqfrom[$i] > $infernal::hitsqto[$i])
	{
	    #hit to reverse strand of query
	    $orient = 1;
	}
	else
	{
	    $orient = 0;
	}
	if($use_evalues)
	{
	    printf("%-24s %-6f %d %d %d\n", $infernal::targname_byhit[$i], $infernal::hitevalue[$i], $infernal::hitsqfrom[$i], $infernal::hitsqto[$i], $orient); 
	}
	else
	{
	    printf("%-24s %-6f %d %d %d\n", $infernal::targname_byhit[$i], $infernal::hitbitscore[$i], $infernal::hitsqfrom[$i], $infernal::hitsqto[$i], $orient); 
	}
    }
    
}
1;
