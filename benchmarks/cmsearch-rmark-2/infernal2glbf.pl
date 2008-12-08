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
#    -E <x>         : use E values [default], sets max E-val to keep as <x> [default = 10]
#    -B <x>         : use bit scores, sets min score to keep as <x>
#    -S             : sort scores
#
# SRE, Wed Oct 28 14:03:52 1998
# RCS $Id$
#

use Getopt::Std;
use infernal;

$use_evalues   = 1;
$use_bitscores = 0;
$e_cutoff =   10;
$b_cutoff = 0.0;
$sort_scores = 0;

getopts('E:B:S');
if (defined $opt_E) { $e_cutoff = $opt_E; }
if (defined $opt_B) { $b_cutoff = $opt_B; $use_evalues = 0; $use_bitscores = 1; }
if (defined $opt_S) { $sort_scores = 1; }

$output = join("",<>);
&infernal::ParseINFERNAL($output);

# First determine if infernal was run with or without E-values
$at_least_one_hit = 0;
for($c = 0; $c < $infernal::ncm; $c++)
{
    if($infernal::nhit[$c] > 0)
    {
	$at_least_one_hit = 1;
	if (exists($infernal::hitevalue[$c][0]))
	{
	    $has_evalues = 1;
	}
	else
	{
	    $has_evalues = 0;
	}
	last;
    }
}
if(!($at_least_one_hit))
{
    die("No hits found. Exiting.\n");
}
if($use_evalues && (!$has_evalues))
{
    die("ERROR, trying to use E-values but none reported.\n");
}

# if we're not sorting the scores, print them out in the order cmsearch reported them 
if(!($sort_scores))
{
    for ($c = 0; $c < $infernal::ncm; $c++)
    {
	for ($i = 0; $i < $infernal::nhit[$c]; $i++)
	{
	    if ((($use_bitscores) && $infernal::hitbitscore[$c][$i] > $b_cutoff) ||
		(($use_evalues)  && $infernal::hitevalue[$c][$i]   < $e_cutoff))
	    {
		#printf("%-24s %-6f\n", $infernal::targname[$i], $infernal::seqbitscore{$infernal::targname[$i]}); 
		if($infernal::hitsqfrom[$c][$i] > $infernal::hitsqto[$c][$i])
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
		    printf("%-24s %g %d %d %d\n", $infernal::targname_byhit[$c][$i], $infernal::hitevalue[$c][$i], $infernal::hitsqfrom[$c][$i], $infernal::hitsqto[$c][$i], $orient); 
		}
		else
		{
		    printf("%-24s %-6f %d %d %d\n", $infernal::targname_byhit[$c][$i], $infernal::hitbitscore[$c][$i], $infernal::hitsqfrom[$c][$i], $infernal::hitsqto[$c][$i], $orient); 
		}
	    }
	}
    }
}
else # sort scores 
{
    for ($c = 0; $c < $infernal::ncm; $c++)
    {
	for ($i = 0; $i < $infernal::nhit; $i++)
	{
	    $key = "$c:$i";
	    if($use_evalues)
	    {
		$sc_H{$key} = $infernal::hitevalue[$c][$i];
	    }
	    else
	    {
		$sc_H{$key} = $infernal::hitbitscore[$c][$i];
	    }
	}
	if($use_evalues)
	{
	    @sorted_i_A = sort { $sc_H{$a} <=> $sc_H{$b} } (keys (%sc_H));
	}
	else
	{
	    @sorted_i_A = sort { $sc_H{$b} <=> $sc_H{$a} } (keys (%sc_H));
	}
	for ($j = 0; $j < scalar(@sorted_i_A); $j++)
	{
	    $key = $sorted_i_A[$j];
	    ($c, $i) = split(":", $key);
	    if ((($use_bitscores) && $infernal::hitbitscore[$c][$i] > $b_cutoff) ||
		(($use_evalues)  && $infernal::hitevalue[$c][$i]   < $e_cutoff))
	    {
		#printf("%-24s %-6f\n", $infernal::targname[$c][$i], $infernal::seqbitscore{$infernal::targname[$c][$i]}); 
		if($infernal::hitsqfrom[$c][$i] > $infernal::hitsqto[$c][$i])
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
		    printf("%-24s %g %d %d %d\n", $infernal::targname_byhit[$c][$i], $infernal::hitevalue[$c][$i], $infernal::hitsqfrom[$c][$i], $infernal::hitsqto[$c][$i], $orient); 
		}
		else
		{
		    printf("%-24s %-6f %d %d %d\n", $infernal::targname_byhit[$c][$i], $infernal::hitbitscore[$c][$i], $infernal::hitsqfrom[$c][$i], $infernal::hitsqto[$c][$i], $orient); 
		}
	    }
	}
    }
}
1;


