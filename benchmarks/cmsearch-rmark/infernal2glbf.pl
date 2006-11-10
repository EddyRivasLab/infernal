#! /usr/bin/perl

# @LICENSE@
# EPN 09.06.05
# modified from hmmer2gdf from SRE

# Usage: perl infernal2glbf.pl <cmsearch output file>
#
# Converts cmsearch output to GLBF.
# GLBF format is just like GLF format but with bounds of hits
# <seq name> <score> <start posn> <end posn> <orientation 0 (forward) of 1 (reverse)>
# Order is sorted by bit score, just like original output.
#
# Options from :
#    -E <x>         : sets E value cutoff (globE) to <x>
#    -T <x>         : sets T bit threshold (globT) to <x>
#
# SRE, Wed Oct 28 14:03:52 1998
# RCS $Id$
#

use Getopt::Std;
use infernal;

$globE =  999999.;
$globT = 0.0;

getopts('E:T:');
if (defined $opt_E) { $globE     = $opt_E; }
if (defined $opt_T) { $globT     = $opt_T; }

$output = join("",<>);
&infernal::ParseINFERNAL($output);

for ($i = 0; $i < $infernal::nhit; $i++)
{
    if ($infernal::hitbitscore[$i] > $globT)
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
	printf("%-24s %-6f %d %d %d\n", $infernal::targname_byhit[$i], $infernal::hitbitscore[$i], $infernal::hitsqfrom[$i], $infernal::hitsqto[$i], $orient); 
    }
}


