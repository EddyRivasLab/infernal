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

# if we're not sorting the scores, print them out in the order cmsearch reported them 
if(!($sort_scores))
{
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
}
else # sort scores 
{
    for ($i = 0; $i < $infernal::nhit; $i++)
    {
	if($use_evalues)
	{
	    $sc_H{$i} = $infernal::hitevalue[$i];
	}
	else
	{
	    $sc_H{$i} = $infernal::hitbitscore[$i];
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
	$i = $sorted_i_A[$j];
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
}
1;


