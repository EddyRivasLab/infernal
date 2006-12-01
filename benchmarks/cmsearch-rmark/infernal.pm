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
#    printf "The top scoring sequence is %s\n", $infernal:targname[0];
#    printf "The total number of hits is %d\n", $infernal::nhit;
#
# Data made available:
#    $ntarget      - total number of sequences with at least 1 hit
#    @targname     - array of target names
#    %seqnhit      - number of sequences hit (indexed by target name)
#
#    $nhit         - total number of hits
#    @hitname      - target names hit
#    @hitnum       - hit number (starts at 0)
#    @hitsqfrom    - sequence from coords (start positions)
#    @hitsqto      - sequence to coords (end positions)
#    @hitbitscore  - hit bit scores
#    @hitevalue    - hit E-values (only available if --E option used)
#    @hitpvalue    - hit P-values (only available if --E option used)

# Data not made available that a future cmsearch might output:
#
#    $query        - name of query CM or sequence
#    $querydesc    - description of query,
#
#    %targdesc     - target descriptions (indexed by target name)
#    %seqscore     - per-seq score or bit score (indexed by target name)
#    %seqeval      - per-seq E-value (indexed by target name) NOT YET IMPLEMENTED IN INFERNAL
#
#    @hitsqbounds  - e.g. "[]" or ".." for seq
#    @hitcmto     - array of cm-to coords
#    @hitcmfrom   - array of cm-from coords
#    @hitcmbounds - e.g. "[]" or ".." for CM
#
#    $aligndata    - the raw alignment text (currently not parsed further)
#
sub ParseINFERNAL {
    my($output) = @_;
    #my($inhit, $inseq, $inali);
    my(@lines, $line);

    #$query       = "";
    #$querydesc   = "";

    $ntarget     = 0;
    @targname    = ();
    @targname_byhit = ();
    #targdesc    = ();
    #%seqscore    = ();
    #%seqeval     = ();
    %seqnhit     = ();

    $nhit        = 0;
    @hitname     = ();
    @hitnum      = ();
    @hitsqfrom   = ();
    @hitsqto     = ();
    #@hitsqbounds = ();
    #@hithmmfrom  = ();
    #@hithmmto    = ();
    #@hithmmbounds= ();
    @hitbitscore    = ();
    @hitevalue   = ();
    @hitpvalue   = ();
    #$aligndata   = "";

    @lines = split(/^/, $output);
    $nhit=0;
    $ntarget=0;
    foreach $line (@lines) 
    {
	chomp $line;
	if ($line =~ /^sequence:\s+(.+)/)
	{
	    $targname[$ntarget] = $1;
	    $ntarget++;
	}
	# if statistics (E and P values) not reported:
	elsif ($line =~ /^hit\s+(\d+)\s+\:\s+(\d+)\s+(\d+)\s+(\S+)\s+bits\S*$/)
	{
	    $hitnum[$nhit]      = $1;
	    $hitsqfrom[$nhit]   = $2; 
	    $hitsqto[$nhit]     = $3;
	    $hitbitscore[$nhit] = $4;
	    $targname_byhit[$nhit]=$targname[$ntarget-1];
	    $nhit++;
	}
        elsif ($line =~ /^hit\s+(\d+)\s+\:\s+(\d+)\s+(\d+)\s+(\S+)\s+bits\s+E\s+\=\s+(\S+)\,\s+P\s+\=\s+(\S+)\s*$/)
	{
	    $hitnum[$nhit]      = $1;
	    $hitsqfrom[$nhit]   = $2; 
	    $hitsqto[$nhit]     = $3;
	    $hitbitscore[$nhit] = $4;
	    $hitevalue[$nhit]   = $5;
	    $hitpvalue[$nhit]   = $6;
	    $targname_byhit[$nhit]=$targname[$ntarget-1];
	    $nhit++;
	}
	1;
    }
}
1;
