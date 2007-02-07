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

use constant FASTA_LINE_LENGTH       => 50;
use Getopt::Std;
use infernal;

$use_evalues   = 1;
$use_bitscores = 0;
$e_cutoff =   10;
$b_cutoff = 0.0;
$sort_scores = 1;
$do_extract = 0;

getopts('E:B:S:X:');
if (defined $opt_X) { $do_extract = 1; $db_file = $opt_X; }
if (defined $opt_E) { $e_cutoff = $opt_E; }
if (defined $opt_B) { $b_cutoff = $opt_B; $use_evalues = 0; $use_bitscores = 1; }
if (defined $opt_S) { $sort_scores = 1; }

$usage = "Usage: perl cmsearch_pp.pl\n\t<cmsearch output file>\n\t<id for output files>\n";
$options_usage  = "\nOptions:\n\t";
$options_usage .= "-X <f> : extract hit subseqs from database in <f>\n\t";
$options_usage .= "-E <x> : use E values [default], sets max E-val to keep as <x> [df= 10]\n\t";
$options_usage .= "-B <x> : use bit scores, sets min score to keep as <x>\n";

if(scalar(@ARGV) != 2)
{
    #printf("argv: %d\n", scalar(@ARGV));
    print $usage;
    print $options_usage;
    exit();
}

($cmsearch_output, $id) = @ARGV;

 open(IN, $cmsearch_output);
$output = join("",<IN>);
 close(IN);
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

# Get all scores
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

# Sort by score and print out info about hits above threshold
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
	$gc_content = $infernal::hitgccontent[$i];
	#printf("%-24s %-6f\n", $infernal::targname[$i], $infernal::seqbitscore{$infernal::targname[$i]}); 
	if($infernal::hitsqfrom[$i] > $infernal::hitsqto[$i])
	{
	    #hit to reverse strand of query
	    $orient = 1;
	    $start = $infernal::hitsqto[$i] - 1; 
	    $end   = $infernal::hitsqfrom[$i] - 1; 
	    # cmsearch reports first res as posn 1, but it's posn 0 in the string
	}
	else
	{
	    $orient = 0;
	    $start = $infernal::hitsqfrom[$i] - 1; 
	    $end   = $infernal::hitsqto[$i] - 1; 
	    # cmsearch reports first res as posn 1, but it's posn 0 in the string
	}
	
	$targname = $infernal::targname_byhit[$i];
	#printf("targname: $targname orient: $orient start: $start\n");
	push(@{$targ_starts_AHA[$orient]{$targname}}, $start);
	$targ_ends_AHH[$orient]{$targname}{$start} = $end;
	$targ_gc_AHH[$orient]{$targname}{$start} = $gc_content;
	
	if($use_evalues)
	{
	    $score = "E" . $infernal::hitevalue[$i];
	    $targ_scores_AHH[$orient]{$targname}{$start} = $score;
	    printf("%-24s %-6f %d %d %d (GC=%d)\n", $infernal::targname_byhit[$i], $infernal::hitevalue[$i], $infernal::hitsqfrom[$i], $infernal::hitsqto[$i], $orient, $gc_content); 
	}
	else
	{
	    $score = "B" . $infernal::hitbitscore[$i];
	    $targ_scores_AHH[$orient]{$targname}{$start} = $score;
	    printf("%-24s %-6f %d %d %d (GC=%d)\n", $infernal::targname_byhit[$i], $infernal::hitbitscore[$i], $infernal::hitsqfrom[$i], $infernal::hitsqto[$i], $orient, $gc_content); 
	}
	$subseq_name = $targname . "|" . ($start+1) . "-" . ($end+1) . "|" . $orient . "|" . $gc_content . "|" . $score;
	push(@subseq_rank_A, $subseq_name);
    }
}


############################################
# Optionally, extract the hits from the db # 
############################################
if($do_extract)
{
    for($t = 0; $t < $infernal::ntarget; $t++)
    {
	for($orient = 0 ; $orient <= 1; $orient++)
	{
	    # targ_starts_AHA[]{$key}[] = is an array of starting points of hits
	    #                            in target seq $key
	    # targ_ends_AHH[]{$key1}{$key2} = is the end point of hit that starts
	    #    
	    # targ_scores_AHH[]{$key1}{$key2} = is bit score or E value of the hit that
	    #                                starts at $key2 in target seq $key1
	    $targ = $infernal::targname[$t];
	    @sorted_starts = sort {$a <=> $b} @{$targ_starts_AHA[$orient]{$targ}};
	    foreach $start (@sorted_starts)
	    {
		#printf("$targ $start %d %s\n", $targ_ends_AHH[$orient]{$targ}{$start}, 
		 #      $targ_scores_AHH[$orient]{$targ}{$start});
	    }
	}
    }
    for($orient = 0; $orient <= 1; $orient++)
    {
	read_fasta_saving_subseqs_efficiently($db_file, \%{$subseq_H[$orient]}, 
					      \%{$targ_starts_AHA[$orient]}, \%{$targ_ends_AHH[$orient]}, 
					      \%{$targ_gc_AHH[$orient]},
					      \%{$targ_scores_AHH[$orient]}, $orient);
    }

    # print out the sequences
    $subseq_file = $id . "_hits.fa";
    for($orient = 0; $orient <= 1; $orient++)
    {
	foreach $subseq (keys(%{$subseq_H[$orient]}))
	{
	    #printf("subseq: $subseq\n");
	    $subseq_all_H{$subseq} = $subseq_H[$orient]{$subseq};
	}
    }
    
    print_fasta_ordered(\%subseq_all_H, $subseq_file, \@subseq_rank_A);
    $num_hits = scalar(@subseq_rank_A);
    printf("*********************************************************\n");
    printf("* Output file: $subseq_file created with top -%5d    \n", $num_hits);
    printf("*              hits from %-20s           \n", ($cmsearch_output . "."));
    printf("*********************************************************\n");
}
1;

#################################################################
# subroutine : read_fasta_saving_subseqs_efficiently
# sub class  : crw and sequence
# 
# EPN, Sun Dec 17 19:07:58 2006
#
# purpose : Open and read a fasta file, saving subsequences to
#           a sequence hash as dictacted by input information.
#
# args : (1) $in_file
#            name of .fa file in current directory
#        (2) $subseq_HR
#            reference to the hash that will contain the subsequence
#            information.  Fasta description line used as key for
#            each sequence, sequence is value.
#        (3) $starts_HAR
#            ref to hash of arrays, key is seq name, array is
#            start positions of subseqs to extract
#        (4) $ends_HHR
#            ref to hash of hashes, 1D key is seq name, 2D
#            key is start point
#        (5) $gc_HHR
#            ref to hash of hashes, 1D key is seq name, 2D
#            key is start point, value is gc content
#        (6) $add2names_HHR
#            ref to hash of hashes, 1D key is seq name, 2D
#            key is string to add to name of subseq
#        (7) $do_revcomp
#            TRUE to reverse complement the sequences, false
#            not to.
################################################################# 
sub read_fasta_saving_subseqs_efficiently
{
    #printf("in read_fasta_saving_subseqs_efficiently\n");
    if(scalar(@_) != 7)
    {
	die("ERROR read_fasta_saving_subseqs_efficiently takes exactly 7 arguments\n");
    }
    ($in_file, $subseq_HR, $starts_HAR, $ends_HHR, $gc_HHR, $add2names_HHR, $do_revcomp) = @_;
    open(IN, $in_file);
    #printf("in read_fasta do_revcomp $do_revcomp\n");
    
    #chomp up beginning blank lines
    $line = <IN>;
    while(!($line =~ m/^>/))
    {
	$line = <IN>;
	print("$line\n");
    }
    chomp $line;
    $targname = $line;
    $targname =~ s/^>//;
    $targname =~ s/\s*$//;
    $targname =~ s/\s+.+$//;
    @sorted_starts = sort {$a <=> $b} @{$starts_HAR->{$targname}};
    
    foreach $start (@sorted_starts)
    {
	#printf("start: $start end: %d\n", $ends_HHR->{$targname}{$start});
    }
    $subseq_idx = 0;
    $seq_posn = 0;
    if(scalar(@sorted_starts) > 0)
    {
	$subseq_remains = 1;
	$next_start = $sorted_starts[$subseq_idx];
	$next_end   = $ends_HHR->{$targname}{$next_start};
	$subseq_name = $targname . "|" . ($next_start+1) . "-" . ($next_end+1) . "|" . $do_revcomp . "|" . $gc_HHR->{$targname}{$next_start} . "|" . $add2names_HHR->{$targname}{$next_start};
    }

    $line_ct = 0;
    while($line = <IN>)
    {
	$line_ct++;
	chomp $line;
	while((!($line =~ m/^>/)) && ($line ne ""))
	{
	    if(!($donotread))
	    {
		$seq_posn += length($line);
	    }
	    #printf("seq_posn: $seq_posn\n");
	    
	    if($subseq_remains)
	    {
		if(!($in_subseq) && $seq_posn < $next_start)
		{
		    #printf("case1: doing nothing\n");
		    $donotread = 0;
		    ;#do nothing
		}
		elsif(!($in_subseq) && $seq_posn >= $next_start)
		{
		    #our subseq starts (and may end) in this chunk
		    #printf("case2: subseq starts here\n");
		    $prev_seq_posn = $seq_posn - length($line);
		    #printf("prev seq posn: $prev_seq_posn next_start: $next_start\n");
		    $to_add = $line;
		    if($next_start != $prev_seq_posn)
		    {
			$to_add = substr($to_add, ($next_start - $prev_seq_posn));
		    }
		    #printf("to_add: $to_add\n");
		    $in_subseq = 1;
		    $donotread = 0;
		    
		    if(($next_end+1) <= $seq_posn)
		    {
			#printf("case 2 seq ends here also\n");
			if(($next_end + 1) != $seq_posn)
			{
			    #printf("next_end:%d - seq_posn:%d: %d\n", $next_end, $seq_posn, ($next_end-$seq_posn));
			    $to_add = substr($to_add, 0, ($next_end-$seq_posn+1));
			}
			#reset our flags
			$subseq_idx++;
			if($subseq_idx < scalar(@sorted_starts))
			{
			    $subseq_HR->{$subseq_name} .= $to_add;
			    $next_start = $sorted_starts[$subseq_idx];
			    $next_end   = $ends_HHR->{$targname}{$next_start};
			    $in_subseq = 0;
			    $subseq_name = $targname . "|" . ($next_start+1) . "-" . ($next_end+1) . "|" . $do_revcomp . "|" . $gc_HHR->{$targname}{$next_start} . "|" . $add2names_HHR->{$targname}{$next_start};
			    # the next subseq may start in this line, so we flag not to read another line yet
			    $donotread = 1;
			    #printf("new subseq_idx: $subseq_idx start: $next_start end: $next_end $subseq_name\n");
			}
			else
			{
			    $subseq_HR->{$subseq_name} .= $to_add;
			    $subseq_remains = 0;
			    $in_subseq = 0;
			}
		    }
		    else
		    {
			$subseq_HR->{$subseq_name} .= $to_add;
		    }			
		}
		elsif($in_subseq && $seq_posn <= $next_end)
		{
		    #our subseq includes this full chunk
		    $subseq_HR->{$subseq_name} .= $line;
		    $donotread = 0;
		}
		elsif($in_subseq && $seq_posn > $next_end)
		{
		    #our subseq ends in this chunk
		    #printf("subseq ends here seq_posn: $seq_posn, next_end: $next_end\n");
		    $to_add = $line;
		    if(($next_end + 1) != $seq_posn)
		    {
			$to_add = substr($to_add, 0, ($next_end-$seq_posn+1));
		    }
		    $subseq_HR->{$subseq_name} .= $to_add;
		    
		    #reset our flags
		    $subseq_idx++;
		    if($subseq_idx < scalar(@sorted_starts))
		    {
			$next_start = $sorted_starts[$subseq_idx];
			$next_end   = $ends_HHR->{$targname}{$next_start};
			$in_subseq = 0;
			$subseq_name = $targname . "|" . ($next_start+1) . "-" . ($next_end+1) . "|" . $do_revcomp . "|" . $gc_HHR->{$targname}{$next_start} . "|" . $add2names_HHR->{$targname}{$next_start};
			# the next subseq may start in this line, so we flag not to read another line yet
			$donotread = 1;
			#printf("new subseq_idx: $subseq_idx start: $next_start end: $next_end $subseq_name\n");
		    }
		    else
		    {
			$subseq_remains = 0;
			$in_subseq = 0;
			$donotread = 0;
		    }
		}	
		#printf("curr seq $subseq_name: $subseq_HR->{$subseq_name}\n");
	    }
	    if(!($donotread)) 
	    {
		#printf("donotread FALSE, reading line\n");
		$line = <IN>;
		$line_ct++;
		chomp $line;
	    }
	}
	chomp $line;
	$targname = $line;
	$targname =~ s/^>//;
	$targname =~ s/\s*$//;
	@sorted_starts = sort {$a <=> $b} @{$starts_HAR->{$targname}};
	
	#printf("\nnew targ seq: $targname\n");
	foreach $start (@sorted_starts)
	{
	    #printf("start: $start end: %d\n", $ends_HHR->{$targname}{$start});
	}
	$subseq_idx = 0;
	$seq_posn = 0;
	if(scalar(@sorted_starts) > 0)
	{
	    $subseq_remains = 1;
	    $next_start = $sorted_starts[$subseq_idx];
	    $next_end   = $ends_HHR->{$targname}{$next_start};
	    $subseq_name = $targname . "|" . ($next_start+1) . "-" . ($next_end+1) . "|" . $do_revcomp . "|" . $gc_HHR->{$targname}{$next_start} . "|" . $add2names_HHR->{$targname}{$next_start};
	}
    }
    # reverse complement if necessary and do some checks 
    foreach $subseq (keys(%{$subseq_HR}))
    {
	if($do_revcomp == 1)
	{
	    $revcomp = rev_comp($subseq_HR->{$subseq});
	    $subseq_HR->{$subseq} = $revcomp
	}

	$length = length($subseq_HR->{$subseq});
	$temp = $subseq;
	$temp =~ /.+\|(\d+)\-(\d+)\|[01]\|.+$/;
	$start = $1;
	$end   = $2;
	#printf("subseq: $subseq start: $start end: $end\n");
	if($length != ($end-$start+1))
	{
	    printf("ERROR subseq: $subseq length should be %d, but it's $length\n", ($end-$start+1));
	    die;
	}
	else
	{
	    #printf("checked $subseq length: $length\n");
	}
    }
    #debug_print_hash($subseq_HR, "subseq");
    #trim_keys_in_hash($seq_hash_ref, DEFAULT_MAX_SEQ_HEADER_LENGTH);
}

#################################################################
# subroutine : rev_comp
# sub class  : sequence
#
# EPN 09.13.05
#
# purpose : Given a string that is a RNA sequence string, return the
#           reverse complement.  String may have gaps.
#
# args : (1) $seq
#            seq string to rev comp
################################################################# 
sub rev_comp
{
    ($seq) = $_[0];

    #print("in rc in seq is :\n$seq\n");
    $comp = $seq;
    #$comp =~ tr/AUCGaucg/UAGCuagc/;
    $comp =~ tr/ACGUTRYMKSWHDBVacgutrymkswhdbv/UGCAAYRKMSWDHVBugcaayrkmswdhvb/;
    #leave X's, N's and gaps
    @comp_arr = split("", $comp);
    @rev_comp_arr = reverse @comp_arr;
    $rc = "";
    foreach $letter (@rev_comp_arr)
    {
	$rc .= $letter;
    }
    #print("in rc returning rc is :\n$rc\n");
    return $rc;
}
#################################################################

#################################################################
# subroutine : print_fasta_ordered
# sub class  : sequence
# 
# EPN 04.20.05
# 
# purpose : Create a single new file and print all sequences
#           in the input hash to that file in the order 
#           specified in the array referenced in $aln_order_arr_ref
#
# args : (1) $seq_hash_ref 
#            reference to hash with sequences as values, and
#            headers as keys
#        (2) $out_file_name
#            name of resulting gapless fasta file
#        (3) $aln_order_arr_ref
#            reference to an array with order of seq headers
#            we want alignment to be printed in
################################################################# 

sub print_fasta_ordered
{
    ($seq_hash_ref, $out_file_name, $aln_order_arr_ref) = @_;
    
    #open and close file just so to ensure its empty
    open(OUT, ">" . $out_file_name);
    close(OUT);
    
    foreach $header (@{$aln_order_arr_ref})
    {
	if(!(exists($seq_hash_ref->{$header})))
	{
	    die "ERROR in print_fasta_ordered - no sequence with key $header that exists in the alignment order array\n";
	}
	$curr_seq = $seq_hash_ref->{$header};
	print_fasta_single_seq($out_file_name, ">>", $header, $curr_seq);
    }
}
#################################################################
# subroutine : print_fasta_single_seq
# sub class  : sequence
# 
# EPN 03.08.05
# 
# purpose : Print a single fasta sequence to a specified file
#           using line lengths of FASTA_LINE_LENGTH (a constant)
#
# args : (1) $out_file_name
#            name of file to open, print to, and close
#        (2) $output_mode
#            either ">" or ">>"
#        (3) $header
#            sequence header to be placed after ">"
#        (4) $seq
#            sequence to print (after breaking up into multiple lines)
################################################################# 

sub print_fasta_single_seq
{
    ($out_file_name, $output_mode, $header, $seq) = @_;
    #print("in print fasta seq\n");
    #print("out file name is $out_file_name\n");
    #print("output mode is $output_mode\n");
    #print("header is $header\n");
    #print("seq is $seq\n");
    
    #open file with output mode
    open(OUT, $output_mode . $out_file_name);
    if(length($seq) != 0)
    {
	
	print OUT (">" . $header . "\n");
	
	$index = 0;
	$length = length($seq);
	while($index < $length)
	{
	    $substr = substr($seq, $index, FASTA_LINE_LENGTH);
	    print OUT ($substr . "\n");
	    $index += FASTA_LINE_LENGTH;
	}
	close(OUT);
    }
}
