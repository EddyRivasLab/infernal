#! /usr/bin/perl

# cmsearch_pp.pl
# EPN, Sun Dec 17 18:01:06 2006
#
# Usage: cmsearch_pp.pl [options] <cmsearch output file> <id for output files>
#
# Parses cmsearch output and optionally extracts subseqs from the database 
# searched.
#
# Also prints out GLBF format:
# GLBF format is just like GLF format but with bounds of hits
# <seq name> <score> <start posn> <end posn> <orientation 0 (forward) of 1 (reverse)>
# Score is either E-value (if cmsearch was run with stats enabled), else it's bit score 
# Order is sorted by sequence position. 
#
# Options from :
#    -X <f>         : extract hit subseqs from database in <f>
#    -E <x>         : use E values [default], sets max E-val to keep as <x> [df = 10]
#    -B <x>         : use bit scores, sets min score to keep as <x>
#    -A <x>         : sort all hits across CMs [default: sort hits for each CM]
#    -Q <x>         : print top <x> hits per query
#    -T <x>         : print top <x> hits per target

require "sre.pl";
use constant FASTA_LINE_LENGTH       => 50;
use Getopt::Std;
use infernal;

$use_evalues   = 1;
$use_bitscores = 0;
$e_cutoff =   10;
$b_cutoff = 0.0;
$sort_scores = 1;
$sort_all_scores = 0;
$do_extract = 0;
$do_overlap = 0;
$do_top_query = 0;
$do_top_target = 0;

getopts('E:B:X:AR:Q:T:');
if (defined $opt_X) { $do_extract = 1; $db_file = $opt_X; }
if (defined $opt_E) { $e_cutoff = $opt_E; }
if (defined $opt_B) { $b_cutoff = $opt_B; $use_evalues = 0; $use_bitscores = 1; }
if (defined $opt_A) { $sort_all_scores = 1; }
if (defined $opt_R) { $remove_overlaps = 1; $overlap_fraction = $opt_R; }
if (defined $opt_Q) { $do_top_query = 1;  $ntop_query = $opt_Q; }
if (defined $opt_T) { $do_top_target = 1; $ntop_target= $opt_T; }

$usage = "Usage: perl cmsearch_pp.pl\n\t<cmsearch output file>\n\t<output file (ONLY if -X enabled)>\n";
$options_usage  = "\nOptions:\n\t";
$options_usage .= "-X <f> : extract hit subseqs from database in <f> using 'sfetch'\n\t";
$options_usage .= "-E <x> : use E values [default], sets max E-val to keep as <x> [df= 10]\n\t";
$options_usage .= "-B <x> : use bit scores, sets min score to keep as <x>\n\t";
$options_usage .= "-A     : sort all hits across CMs [default: sort hits for each CM]\n\t";
$options_usage .= "-R <x> : remove overlapping hits of <x> fraction (0: no overlap allowed)\n\t";
$options_usage .= "-Q <x> : print top <x> hits per query\n\t";
$options_usage .= "-T <x> : print top <x> hits per target\n\n";
#    -T <x>         : print top <x> hits per target

if(scalar(@ARGV) == 2)
{
    if(!($do_extract))
    {
	print $usage;
	print $options_usage;
	exit();
    }
    ($cmsearch_output, $extract_out) = @ARGV;
}
elsif(scalar(@ARGV) == 1)
{
    if($do_extract)
    {
	print $usage;
	print $options_usage;
	exit();
    }
    ($cmsearch_output) = $ARGV[0];
}
else
{
    print $usage;
    print $options_usage;
    exit();
}
if($do_extract)
{
    if(! (-e "$db_file")) { die("ERROR, -X $db_file enabled but $db_file does not exist.\n") };
}

 open(IN, $cmsearch_output);
 $output = join("",<IN>);
 close(IN);
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
    die("No hits found in $cmsearch_output. Exiting.\n");
}
if($use_evalues && (!$has_evalues))
{
    die("ERROR, trying to use E-values but none reported.\n");
}

# Get all scores
if($use_evalues)
{
    for ($c = 0; $c < $infernal::ncm; $c++)
    {
	for ($i = 0; $i < $infernal::nhit[$c]; $i++)
	{
	    $all_key = "$c:$i";
	    $all_sc_H{$all_key} = $infernal::hitevalue[$c][$i];
	    $sc_AH[$c]{$i}     = $infernal::hitevalue[$c][$i];
	}
	if(!($sort_all_scores))
	{ 
	    @{$sorted_i_AA[$c]} = sort { $sc_AH[$c]{$a} <=> $sc_AH[$c]{$b} } (keys (%{$sc_AH[$c]}));
	}
    }
    if($sort_all_scores)
    {
	@sorted_c_i_A = sort { $all_sc_H{$a} <=> $all_sc_H{$b} } (keys (%all_sc_H));
    }
}
else #use bit scores
{
    for ($c = 0; $c < $infernal::ncm; $c++)
    {
	for ($i = 0; $i < $infernal::nhit[$c]; $i++)
	{
	    $all_key = "$c:$i";
	    $all_sc_H{$all_key} = $infernal::hitbitscore[$c][$i];
	    $sc_AH[$c]{$i}     = $infernal::hitbitscore[$c][$i];
	}
	if(!($sort_all_scores))
	{ 
	    @{$sorted_i_AA[$c]} = sort { $sc_AH[$c]{$b} <=> $sc_AH[$c]{$a} } (keys (%{$sc_AH[$c]}));
	}
    }
    if($sort_all_scores)
    {
	@sorted_ci_A = sort { $all_sc_H{$b} <=> $all_sc_H{$a} } (keys (%all_sc_H));
    }
}
# print out info about hits above threshold
if(!($sort_all_scores)) # build @sorted_ci_A by concatenating sorted list for each CM
{
    @sorted_ci_A = ();
    for ($c = 0; $c < $infernal::ncm; $c++)
    {
	for ($j = 0; $j < scalar(@{$sorted_i_AA[$c]}); $j++)
	{
	    $i = $sorted_i_AA[$c][$j];
	    push(@sorted_ci_A, ($c. ":" . $i));
	}
    }
}
$prev_c = 0;
if($use_evalues) { $sctype = "E"; } 
else { $sctype = "B"; }
for($x = 0; $x < scalar(@sorted_ci_A); $x++)
{
    $print_lines_A[$x] = 1; #this will be changed to 0 if we're removing overlaps
                            #and determine this line to be have a better scoring
                            #overlap

    ($c, $i) = split(":", $sorted_ci_A[$x]);

    if ((($use_bitscores) && $infernal::hitbitscore[$c][$i] > $b_cutoff) ||
	(($use_evalues)   && $infernal::hitevalue[$c][$i]   < $e_cutoff))
    {
	$gc_content = $infernal::hitgccontent[$c][$i];
	#printf("%-24s %-6f\n", $infernal::targname[$i], $infernal::seqbitscore{$infernal::targname[$i]}); 
	if($infernal::hitsqfrom[$c][$i] > $infernal::hitsqto[$c][$i])
	{
	    #hit to reverse strand of query
	    $orient = 1;
	    $start = $infernal::hitsqto[$c][$i]; 
	    $end   = $infernal::hitsqfrom[$c][$i]; 
	}
	else
	{
	    $orient = 0;
	    $start = $infernal::hitsqfrom[$c][$i]; 
	    $end   = $infernal::hitsqto[$c][$i]; 
	}
	$cm       = $infernal::cm[$c];
	$targname = $infernal::targname_byhit[$c][$i];
	if(!(exists($all_targets_H{$targname}))) { $all_targets_H{$targname} = 1; }
	if($use_evalues) { $sc = $infernal::hitevalue[$c][$i]; }
	else             { $sc = $infernal::hitbitscore[$c][$i]; }
	if($remove_overlaps)
	{
	    $start_key = $start . "." . $c; # CM #3 hit starts at 124 will be 124.3
	    push(@{$targ_starts_AH[$orient]{$targname}}, $start_key);
	    $targ_ends_AHH[$orient]{$targname}{$start_key}    = $end;
	    $targ_scores_AHH[$orient]{$targname}{$start_key}  = $sc;
	    $targ_x_AHH[$orient]{$targname}{$start_key}       = $x;
	}
	if($cm =~ m/\w/)
	{
	    $out_lines_A[$x] = sprintf("%-24s %-24s %7.2f %9d %9d %d (GC=%2d)\n", $cm, $targname, $sc, $start, $end, $orient, $gc_content); 
	    $cm .= "|";
	    if($do_top_query) { $cm_for_out_lines_A[$x]     = $cm;    }
	    if($do_top_target){ $target_for_out_lines_A[$x] = $targname;}
	}
	else
	{
	    $out_lines_A[$x] = sprintf("%-24s %7.2f %9d %9d %d (GC=%2d)\n", $targname, $sc, $start, $end, $orient, $gc_content); 
	    if($do_top_query) { $cm_for_out_lines_A[$x]     = "ONLYONE";}
	    if($do_top_target){ $target_for_out_lines_A[$x] = $targname;  }
	}
	if($do_extract)
	{
	    $sfetch_lines_A[$x] = $targname . ":" . $start . ":" . $end . ":" . $orient . ":" . 
		$cm . $targname . "|" . $start . "-" . $end . "|" . $orient . "|" . $gc_content . "|" . $sctype. $sc;
	}
    }
}

# remove overlaps
if($remove_overlaps)
{
    @all_targets_A = keys(%all_targets_H);
    $nall_targets = scalar(@all_targets_A);
    for($t = 0; $t < $nall_targets; $t++)
    {
	for($orient = 0 ; $orient <= 1; $orient++)
	{
	    $targname = $all_targets_A[$t];
	    #printf("checking target: $targname orient $orient\n");
	    @sorted_starts = sort {$a <=> $b} @{$targ_starts_AH[$orient]{$targname}};
	    for($s = 0; $s < (scalar(@sorted_starts) -1); $s++)
	    {
		$start = $sorted_starts[$s];
		$real_start = $start;
		$real_start =~ s/\.(\d+)$//;
		$end   = $targ_ends_AHH[$orient]{$targname}{$start};
		$sc    = $targ_scores_AHH[$orient]{$targname}{$start};
		$x     = $targ_x_AHH[$orient]{$targname}{$start};
		
		$ns = $s+1;
		#printf("end: $end, next_start: %d\n", $sorted_starts[$ns]);
		while($sorted_starts[$ns] < ($end+1)) #remember $start and $nstart
		                                      #have "\.<cm num>" appended at end
		{
		    $nstart = $sorted_starts[($ns)];
		    $real_nstart = $nstart;
		    $real_nstart =~ s/\.(\d+)$//;
		    $nend   = $targ_ends_AHH[$orient]{$targname}{$nstart};
		    $nsc    = $targ_scores_AHH[$orient]{$targname}{$nstart};
		    $nx     = $targ_x_AHH[$orient]{$targname}{$nstart};
		    
		    $too_much_overlap = overlap($real_start, $end, $real_nstart, $nend, 
						$overlap_fraction);
		    if($too_much_overlap)
		    {
			if($use_evalues && ($sc < $nsc)) # remove $ns
			{ $print_lines_A[$nx] = 0; }
			elsif($use_evalues)              # remove $s
			{ $print_lines_A[$x]  = 0; last; }     
			elsif($use_bitscores && ($sc >= $nsc))  #remove $ns
			{ $print_lines_A[$nx] = 0; }
			elsif($use_bitscores)                   #remove $s
			{ $print_lines_A[$x]  = 0; last; }
		    }
		    $ns++;
		    if($ns == scalar(@sorted_starts)) { last; }
		}
	    }
	}
    }
}
# print output, from which overlapping hits may have been removed (if $remove_overlaps)
$prev_c = 0;
for($x = 0; $x < scalar(@sorted_ci_A); $x++)
{
    ($c, $i) = split(":", $sorted_ci_A[$x]);
    if((!($sort_all_scores)) && ($prev_c != $c)) { printf("\n"); }
    $prev_c = $c;
    if($print_lines_A[$x])
    {
	$print_flag = 1;
	if($do_top_query)
	{
	    if(($printed_per_cm_H{$cm_for_out_lines_A[$x]}++) >= $ntop_query)
	    {
		$print_flag = 0;
	    }
	}
	if($do_top_target)
	{
	    if(($printed_per_target_H{$target_for_out_lines_A[$x]}++) >= $ntop_target)
	    {
		$print_flag = 0;
	    }
	}
	if($print_flag) { printf $out_lines_A[$x]; }
    }
}

###############################################################################
# Optionally, extract the hits (after possibly removing overlaps) from the db # 
###############################################################################
if($do_extract)
{
    @sys_sfetch_lines_A = ();
    @sys_rm_lines_A = ();
    for($x = 0; $x < scalar(@sorted_ci_A); $x++)
    {
	if($print_lines_A[$x])
	{
	    $sfetch_line = $sfetch_lines_A[$x];
	    open(OUT, ">" . $extract_out);
	    close(OUT);
	    ($seq, $start, $end, $orient, $extra) = split(":", $sfetch_line);
	    
	    $extra =~ s/\|/\\\|/g;
	    $extra =~ s/\-/\\\-/g;
	    $tmp  = &tempname(); # requires sre.pl
	    
	    if($extra ne "")
	    {
		$sfetch_com = "sfetch -f $start -t $end -F fasta -r $extra -d $db_file $seq > $tmp";
	    }
	    else
	    {
		$sfetch_com = "sfetch -f $start -t $end -F fasta -d $db_file $seq > $tmp";
	    }
	    if($orient == 1)
	    {
		#set up revcomp calls
		$tmp2  = &tempname(); # requires sre.pl
		$revcomp_com = "revcomp $tmp > $tmp2";
		$cat_com = "cat $tmp2 >> $extract_out";
		$rm_com  = "rm $tmp2";
		push(@sys_sfetch_lines_A, $sfetch_com);
		push(@sys_sfetch_lines_A, $revcomp_com);
		push(@sys_sfetch_lines_A, $cat_com);
		push(@sys_sfetch_lines_A, $rm_com);
	    }
	    else
	    {
		$cat_com = "cat $tmp >> $extract_out";
		$rm_com  = "rm $tmp";
		push(@sys_sfetch_lines_A, $sfetch_com);
		push(@sys_sfetch_lines_A, $cat_com);
		push(@sys_sfetch_lines_A, $rm_com);
	    }
	}
    }
    #open(OUT, ">" . temp);
    foreach $line (@sys_sfetch_lines_A)
    {
	#print OUT ("$line\n");
	system("$line");
    }
    #close(OUT);
}

#################################################################
# subroutine : overlap
#
# EPN 09.15.05
#
# purpose : Determine if one hit overlaps significantly with a given 
#           region by more than $min_overlap_fract.
#
# args : (1) $begin1
#            begin position of region 1
#        (2) $end1
#            end of region 1
#        (3) $begin2 
#            begin position of region 2  
#        (4) $end2 
#            end position of region 2  
#        (5) $min_overlap_fract
#            to return TRUE (1) the overlap between region 1 and 2
#            must be > $min_overlap_fract * len(min(len(region1), len(region2))
################################################################# 
sub overlap
{
    if(scalar(@_) != 5)
    {
	die "ERROR in overlap, exactly 5 arguments expected.\n";
    }

    my($begin1, $end1, $begin2, $end2, $min_overlap_fract) = @_;
    my $overlap = 0;
    my $overlap_fract = 0;
    #four mutually exclusive possibilities of actual overlap
    if(($begin2 <= $begin1) && ($end2 >= $end1))
    {
	$overlap = $end1 - $begin1 + 1;
    }
    elsif(($begin2 > $begin) && ($end2 >= $end))
    {
	$overlap = $end1 - $begin2 + 1;
    }
    elsif(($begin2 <= $begin1) && ($end2 < $end1))
    {
	$overlap = $end2 - $begin1 + 1;
    }
    elsif(($begin2 > $begin1) && ($end2 < $end1))
    {
	$overlap = $end2 - $begin2 + 1;
    }
    $len1 = $end1-$begin1+1;
    $len2 = $end2-$begin2+1;
    $min_len = $len1;
    if($len2 < $len1) { $min_len = $len2; }
    $overlap_fract = $overlap / $min_len;
    if($overlap_fract > $min_overlap_fract)
    {
	#printf("returning 1, overlap_fract: $overlap_fract\n");
	return 1;
    }
    #printf("returning 0, overlap_fract: $overlap_fract\n");
    return 0;
}
#################################################################












