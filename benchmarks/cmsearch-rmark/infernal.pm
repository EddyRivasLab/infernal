# Eric Nawrocki
# 09.06.05
# infernal.pm
# Perl routines for parsing INFERNAL output
# 
# based on hmmer.pm
# SRE, Wed Oct 28 11:27:17 1998
# RCS $Id$ 

package infernal;

#------------ ParseInfernal ------------
#
# Parse cmsearch (v0.55->v0.71) output into
# arrays that we can use in perl scripts.
#
# Illustrative example:
#    use infernal;
#    $output = `cmsearch foo.cm RF0005`;
#    &infernal::ParseINFERNAL($output);
#    printf "The total number of hits is %d\n", $infernal::nhit;
#
# Data made available for each CM $c:
#    $ncm              - number of CMs used to search
#    $ntarget[$c]      - total number of sequences with at least 1 hit
#    @{$targname[$c]}  - array of target names
#    %{$seqnhit[$c]}   - number of sequences hit (indexed by target name)
#
#    $cm[$c]        - name of query CM 
#    $cmdesc[$c]    - description of query CM
#
#    $nhit[$c]            - total number of hits
#    @{$hitname[$c]}      - target names hit
#    @{$hitnum[$c]}       - hit number (starts at 0)
#    @{$hitsqfrom[$c]}    - sequence from coords (start positions)
#    @{$hitsqto[$c]}      - sequence to coords (end positions)
#    @{$hitbitscore[$c]}  - hit bit scores
#    @{$hitevalue[$c]}    - hit E-values (only available if -E option used)
#    @{$hitpvalue[$c]}    - hit P-values (only available if -E option used)
#    @{$hitgccontent[$c]} - GC content (an integer 0..100) 
#    @{$hitcmfrom[$c]}    - array of cm-from coords
#    @{$hitcmto[$c]}      - array of cm-to coords

# Data not made available that a future cmsearch might output:
#
#    @hitsqbounds  - e.g. "[]" or ".." for seq
#
#    @hitcmbounds - e.g. "[]" or ".." for CM
#
#    $aligndata    - the raw alignment text (currently not parsed further)
#
sub ParseINFERNAL {
    my($output) = @_;
    my(@lines, $line);

    $cm[0]       = "";
    $cmdesc[0]   = "";

    @targname             = ();
    @{$targname[0]}       = ();
    @targname_byhit       = ();
    @{$targname_byhit[0]} = ();
    @seqnhit              = ();
    %{$seqnhit[0]}        = ();

    @nhit              = ();
    $nhit[0]           = 0;
    @ntarget           = ();
    $ntarget[0]        = 0;

    @hitname           = ();
    @{$hitname[0]}     = ();
    @hitnum            = ();
    @{$hitnum[0]}      = ();
    @hitsqfrom         = ();
    @{$hitsqfrom[0]}   = ();
    @hitsqto           = ();
    @{$hitsqto[0]}     = ();
    @hitcmfrom         = ();
    @{$hitcmfrom[0]}   = ();
    @hitcmto           = ();
    @{$hitcmto[0]}     = ();
    @hitbitscore       = ();
    @{$hitbitscore[0]} = ();
    @hitevalue         = ();
    @{$hitevalue[0]}   = ();
    @hitpvalue         = ();
    @{$hitpvalue[0]}   = ();
    @hitgccontent      = ();
    @{$hitgccontent[0]}= ();

    @lines = split(/^/, $output);
    $ncm = 0;
    $seen_cm = 0;
    foreach $line (@lines) 
    {
	chomp $line;
	########################################################################################
	# 03.21.07 New cmsearch output prints CM name (and results for potentially multiple CMs)
	########################################################################################
	if ($line =~ /^CM (\d+):\s+(.+)$/)
	{
	    if(!($seen_cm)) { $seen_cm = 1; } #$ncm stays 0 for first CM 
	    else { $ncm++; }
	    $cm[$ncm] = $2;
	    if($ncm >= 1)
	    {	    
		$cmdesc[$ncm] = "";
		$nhit[$ncm] = 0;
		$ntarget[$ncm] = 0;
		@{$targname[$ncm]}       = ();
		@{$targname_byhit[$ncm]} = ();
		%{$seqnhit[$ncm]}        = ();
		@{$hitname[$ncm]}     = ();
		@{$hitnum[$ncm]}      = ();
		@{$hitsqfrom[$ncm]}   = ();
		@{$hitsqto[$ncm]}     = ();
		@{$hitcmfrom[$ncm]}   = ();
		@{$hitcmto[$ncm]}     = ();
		@{$hitbitscore[$ncm]} = ();
		@{$hitevalue[$ncm]}   = ();
		@{$hitpvalue[$ncm]}   = ();
		@{$hitgccontent[$ncm]}= ();
	    }
	}
	elsif ($line =~ /^CM desc:\s+(.+)$/)
	{
	    $cmdesc[$ncm] = $1;
	}
	#########################################################################
	# 12.08.06 OLD cmsearch output 0.72 and before, handled by next 3 elsif's
	#########################################################################
	elsif ($line =~ /^sequence:\s+(.+)/)
	{
	    $targname[$ncm][$ntarget] = $1;
	    $ntarget[$ncm]++;
	}
	# if statistics (E and P values) not reported:
	elsif ($line =~ /^hit\s+(\d+)\s*\:\s+(\d+)\s+(\d+)\s+(\S+)\s+bits\S*$/)
	{
	    $hitnum[$ncm][($nhit[$ncm])]      = $1;
	    $hitsqfrom[$ncm][($nhit[$ncm])]   = $2; 

	    $hitsqto[$ncm][($nhit[$ncm])]     = $3;
	    $hitbitscore[$ncm][($nhit[$ncm])] = $4;
	    $targname_byhit[$ncm][($nhit[$ncm])]=$targname[$ncm][$ntarget-1];
	    $nhit++;
	}
        elsif ($line =~ /^hit\s+(\d+)\s*\:\s+(\d+)\s+(\d+)\s+(\S+)\s+bits\s+E\s+\=\s+(\S+)\,\s+P\s+\=\s+(\S+)\s*$/)
	{
	    $hitnum[$ncm][($nhit[$ncm])]      = $1;
	    $hitsqfrom[$ncm][($nhit[$ncm])]   = $2; 
	    $hitsqto[$ncm][($nhit[$ncm])]     = $3;
	    $hitbitscore[$ncm][($nhit[$ncm])] = $4;
	    $hitevalue[$ncm][($nhit[$ncm])]   = $5;
	    $hitpvalue[$ncm][($nhit[$ncm])]   = $6;
	    $targname_byhit[$ncm][($nhit[$ncm])]=$targname[$ncm][$ntarget-1];
	    $nhit[$ncm]++;
	}
	#########################################################################
	# 12.08.06 New cmsearch output (from RSEARCH) picked up by next 4 elsif's
	#########################################################################
	elsif($line =~ /^>(.+)$/)
	{
	    $targname[$ncm][$ntarget] = $1;
	    $ntarget[$ncm]++;
	}
	elsif($line =~ /^\s+Query\s+\=\s+(\d+)\s+\-\s+(\d+)\,\s+Target\s+\=\s+(\d+)\s+\-\s+(\d+)\s*$/)
	{
	    $hitcmfrom[$ncm][($nhit[$ncm])]   = $1;
	    $hitcmto[$ncm][($nhit[$ncm])]     = $2;
	    $hitsqfrom[$ncm][($nhit[$ncm])]   = $3; 
	    $hitsqto[$ncm][($nhit[$ncm])]     = $4;
	    $targname_byhit[$ncm][($nhit[$ncm])]=$targname[$ncm][$ntarget-1];
	}
	# ^Query line always followed by ^Score line
	# ^Score line either has E and P values or doesn't
	elsif ($line =~ /^\s+Score\s+\=\s+(\S+),\s+GC\s+\=\s+(\d+)\s*$/)
	{
	    # no E or P values reported 
	    $hitbitscore[$ncm][($nhit[$ncm])] = $1;
	    $hitgccontent[$ncm][($nhit[$ncm])]= $2;
	    $nhit[$ncm]++;
	}
	elsif ($line =~ /^\s+Score\s+\=\s+(\S+),\s+E\s+\=\s+(\S+)\,\s+P\s+\=\s+(\S+),\s+GC\s+\=\s+(\d+)\s*$/)
	{
	    # E and P values reported 
	    $hitbitscore[$ncm][($nhit[$ncm])] = $1;
	    $hitevalue[$ncm][($nhit[$ncm])]   = $2;
	    $hitpvalue[$ncm][($nhit[$ncm])]   = $3;
	    $hitgccontent[$ncm][($nhit[$ncm])]= $4;
	    $nhit[$ncm]++;
	}

	####################################################
	# Special section for parsing RSEARCH output,
	# only difference with infernal output is GC content
	# not reported.
	####################################################
	# ^Query line always followed by ^Score line
	# ^Score line either has E and P values or doesn't
	elsif ($line =~ /^\s+Score\s+\=\s+(\S+)\s*$/)
	{
	    # no E or P values reported 
	    $hitbitscore[$ncm][($nhit[$ncm])] = $1;
	    $nhit[$ncm]++;
	}
	elsif ($line =~ /^\s+Score\s+\=\s+(\S+),\s+E\s+\=\s+(\S+)\,\s+P\s+\=\s+(\S+)\s*$/)
	{
	    # E and P values reported 
	    $hitbitscore[$ncm][($nhit[$ncm])] = $1;
	    $hitevalue[$ncm][($nhit[$ncm])]   = $2;
	    $hitpvalue[$ncm][($nhit[$ncm])]   = $3;
	    $nhit[$ncm]++;
	}
	1;
    }
    $ncm++; # account for off-by-one with array indexing
}
1;
