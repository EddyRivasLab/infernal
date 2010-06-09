#!/usr/bin/perl
#
# Eric Nawrocki
# Fri Feb 29 13:30:00 2008
# 
# hmm_generate_only.pl - This script reads in a parameter file which describes an HMM
#                        and it generates a specified number of sequences of a specified
#                        length using the HMM.
#
# Based on hmm_generate_embed.pl, but simpler, there's no embedding.
#           
# Modified from Algorithms for Computational Biology Spring 04 HW4 part 1
#############################################################################
use strict;
use constant FASTA_LINE_LENGTH       => 50;
use constant SPECIAL_LINE_LENGTH     => 75; 
use constant DEFAULT_LINE_LENGTH     => 75;   

my $usage ="perl hmm_generate_only.pl\n\t<params file>\n\t<num output seqs>\n\t<length of output seqs>\n\t<output root>\n\t<seed for RNG>\n";

(@ARGV == 5) || die $usage;

my($params_file, $num_out_seqs, $output_len, $out_root, $seed) = @ARGV;

my %init_p_hash;    #hash for initial state probabilities, key is name of state
                    #value is probability of starting in that state
my %trans_p_hash;   #2D hash for transition probabilities, first key is name of
                    #state i, second key is name of state j, value is probability
                    #of transition from state i to state j
my %emit_p_hash;    #2D hash for emission probabilities, first key is name of 
                    #state i, second key is a letter A from the output alphabet,
                    #value is probability of emitting A from state i
my @st_names_arr;   #array of state names in order they appear in input file
my @out_alph_arr;   #array of characters in output alphabet 
my @output_arr;     #array which holds the output


#read in the parameters from the file
read_params($params_file, \%init_p_hash, \%trans_p_hash, \%emit_p_hash, \@st_names_arr, \@out_alph_arr);

#produce output based on HMM
#produce_output($output_len, \%init_p_hash, \%trans_p_hash, \%emit_p_hash, \@st_names_arr, \@out_alph_arr, \@output_arr, $num_out_seqs, $out_root);

my $out_fa_file = $out_root . ".fa";

#produce output from HMM
generate($output_len, \%init_p_hash, \%trans_p_hash, \%emit_p_hash, \@st_names_arr, \@out_alph_arr, \@output_arr, $num_out_seqs, $out_root, $out_fa_file);

sub read_params
{
    #read in parameter file
    my($params_file, $init_p_hash_ref, $trans_p_hash_ref, $emit_p_hash_ref, $st_names_arr_ref, $out_alph_arr_ref) = @_;
    my $i;
    my $j; 
    my $ikey;
    my $jkey;
    my $key;
    my @trans_p_arr;
    my @emit_p_arr;

    open(PFILE, $params_file);

    #read in state names
    @{$st_names_arr_ref} = split " ", <PFILE>;
    my $num_states = scalar(@{$st_names_arr_ref});
    #read in init probs
    my @init_p_arr = split " ", <PFILE>;
    my $sum = 0.;
    for($i=0;$i<$num_states;$i++)
    {
	$init_p_hash_ref->{$st_names_arr_ref->[$i]} = $init_p_arr[$i];
	$sum += $init_p_arr[$i];
    }
    if(abs($sum - 1.0) > 0.00001) { printf("ERROR init sum: $sum != 1\n"); exit(1); }
    #read in trans probs
    for($i=0;$i<$num_states;$i++)
    {
	@trans_p_arr = split " ", <PFILE>;
	$sum = 0.;
	for($j=0;$j<$num_states;$j++)
	{
	    $trans_p_hash_ref->{$st_names_arr_ref->[$i]}{$st_names_arr_ref->[$j]} = $trans_p_arr[$j];
	    $sum += $trans_p_arr[$j];
	}
	if(abs($sum - 1.0) > 0.00001) { printf("ERROR trans sum st: $i $sum != 1\n"); exit(1); }
    }
    @{$out_alph_arr_ref} = split " ", <PFILE>;
    my $len_alph = scalar(@{$out_alph_arr_ref});
    for($i=0;$i<$num_states;$i++)
    {
	$sum = 0.;
	@emit_p_arr = split " ", <PFILE>;
	for($j=0;$j<$len_alph;$j++)
	{
	    $emit_p_hash_ref->{$st_names_arr_ref->[$i]}{$out_alph_arr_ref->[$j]} = $emit_p_arr[$j];
	    $sum += $emit_p_arr[$j];
	}
	if(abs($sum - 1.0) > 0.00001) { printf("ERROR emit sum st $i: $sum != 1\n"); exit(1); }
    }
}

sub choose_from_dist
{
    #given a probability distribution (array of numbers <= 1, which sum to 1)
    #this function returns the index of the array chosen based on that 
    #probability distribution

    my ($curr_prob_dist_arr_ref) = @_;
    my $i;
    my $rand_num;
    my $char_index_to_add;
    my $curr_threshold;
    my $num_chars = length(@{$curr_prob_dist_arr_ref});
    #get a random number greater than or equal to 0 and less than 1
    $rand_num = rand();
    $char_index_to_add = 0;
    $curr_threshold = $curr_prob_dist_arr_ref->[$char_index_to_add];
    while($curr_threshold < $rand_num)
    {
	$char_index_to_add++;
	$curr_threshold += $curr_prob_dist_arr_ref->[$char_index_to_add];
    }
    #char_index_to_add now equals the index of the chosen char in out_alph_arr
    return $char_index_to_add;
}

sub generate
{
    #Given parameters of an HMM, generate seqs from it
    #calls choose_from_dist iteratively
    my($output_len, $init_p_hash_ref, $trans_p_hash_ref, $emit_p_hash_ref, 
       $st_names_arr_ref, $out_alph_arr_ref, $output_arr_ref, $num_seqs,
       $name_root, $out_fa_file) = @_;
    my $i;
    my $j;
    my $k;
    my $q;
    my $ikey;
    my $jkey;
    my @curr_prob_dist_arr = ();
    my $char;
    my $state;
    my $name;
    my $curr_len;
    my $curr_chrom;
    my $curr_begin;
    my $curr_end;
    my $other_begin;
    my $other_end;
    my $curr_orient;
    my $conflict;
    my $key;
    my %chrom_begin_hash = ();
    my %chrom_end_hash = ();
    my $chrom_num;

    my $curr_state;


    my @curr_sorted_chrom_arr_of_hash = ();
    my %hash = ();
    my @temp_sorted_keys = ();

    my $offset;
    my %to_print_seq_hash = ();
    my @ordered_names = ();

    srand(($seed * 100000)); #time|$$ for "random" seed

    
    #produce output of specified length based on the HMM parameters
    $curr_state = $st_names_arr_ref->[rand(scalar (keys %{$init_p_hash_ref}))];

    # create prob distros for emissions and transitions
    my %emit_prob_dist_HA = ();
    my %trans_prob_dist_HA = ();
    my $next_state;
    foreach $curr_state (@{$st_names_arr_ref}) { 
	@{$emit_prob_dist_HA{$curr_state}} = ();
	@{$trans_prob_dist_HA{$curr_state}} = ();
	foreach $char (@{$out_alph_arr_ref})
	{
	    push(@{$emit_prob_dist_HA{$curr_state}}, $emit_p_hash_ref->{$curr_state}{$char});
	}
	foreach $next_state (@{$st_names_arr_ref})
	{
	    push(@{$trans_prob_dist_HA{$curr_state}}, $trans_p_hash_ref->{$curr_state}{$next_state});
	}
    }

    #for each sequence 
    for($i=0; $i < $num_seqs; $i++)
    {
	$name = $name_root . "_" . ($i + 1);

	#create the output 
	@{$output_arr_ref} = ();
	for($q=0;$q<$output_len;$q++)
	{
	    #emit letter of output alphabet based on probability dist
	    $output_arr_ref->[$q] = $out_alph_arr_ref->[choose_from_dist(\@{$emit_prob_dist_HA{$curr_state}})];
	    #printf("curr_state: $curr_state res: $output_arr_ref->[$q]\n");
	    #transition to new state based on probability dist
	    $curr_state = $st_names_arr_ref->[choose_from_dist(\@{$trans_prob_dist_HA{$curr_state}})];
	}
	#print output
	#for($i=0;$i<$output_len;$i++)
	#{
	#    if(($i != 0) && ($i % 60 == 0)) { print("\n"); }
	#    #print $output_arr_ref->[$i];
	#}
	#print("pushed name to ordered names\n");
	$to_print_seq_hash{$name} = array_to_scalar($output_arr_ref);
	push(@ordered_names, $name);
    }
    
    print_fasta_ordered(\%to_print_seq_hash, $out_fa_file, \@ordered_names);
    print_out_file_notice($out_fa_file, "New fasta file containing the " . (scalar(@ordered_names)) . " sequences generated by HMM explained in $params_file.");
}


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
    my ($seq_hash_ref, $out_file_name, $aln_order_arr_ref) = @_;
    
    my $header;
    my $curr_seq;

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
    my ($out_file_name, $output_mode, $header, $seq) = @_;
    #print("in print fasta seq\n");
    #print("out file name is $out_file_name\n");
    #print("output mode is $output_mode\n");
    #print("header is $header\n");
    #print("seq is $seq\n");
    
    #open file with output mode
    my $index;
    my $length;
    my $substr;

    open(OUT, $output_mode . $out_file_name);
    if(length($seq) != 0)
    {
	
	print OUT (">" . $header . "\n");
	
	$index = 0;
	$length = length($seq);
	while($index < $length)
	{
	    $substr = substr($seq, $index, 50);
	    print OUT ($substr . "\n");
	    $index += FASTA_LINE_LENGTH;
	}
	close(OUT);
    }
}


#################################################################
# subroutine : print_out_file_notice
# sub class  : general 
#
# EPN 03.03.05
# 
# purpose : Print an output file 'notice' to standard output
#           given the name of the output file and a short message
#           describing that output file
# 
# args : (1) $file_name
#            name of file
#        (2) $description
#            description of file
#################################################################
sub print_out_file_notice
{
    my ($file_name, $description) = @_;

    my $char = "*";
    my $spec_line = create_special_line($char, SPECIAL_LINE_LENGTH);
    print("$spec_line");
    print(" Output file notice\n");
    print(" File name   : $file_name\n");
    
    my $new_desc = add_newlines_to_string($description, length(" description : "),
					  length(" description : "));
    print(" description : $new_desc\n");
    print("$spec_line");
}


#################################################################
# subroutine : add_newlines_to_string
# sub class  : general
#
# EPN 03.03.05
# 
# purpose : Given a string, create a new string that is identical
#           but has new lines inserted so that the line doesn't wrap
#           when printed to the screen
# 
# args : (1) $string
#            string to add new lines to
#        (2) $indent_len
#            number of spaces to put in indenting in all rows but first
#        (3) $first_line_length
#            length to remove from first line's length
#        (4) $char_to_start
#            character to put at beginning of each line 
#            can be ommitted.
#################################################################
sub add_newlines_to_string
{
    my ($string, $indent_len, $first_line_length, $char_to_start) = @_;
    
    my @string_arr = split(" ", $string);
    my $line_beg = 1;
    my $new_string = $char_to_start;
    my $curr_line_len = $first_line_length;
    my $indent = "";
    my $tok_len;
    my $tok;
    my $i;
    for($i = 0; $i < $indent_len; $i++)
    {
	$indent .= " ";
    }

    foreach $tok (@string_arr)
    {
	$tok_len = length($tok);
	if($tok_len > 0)
	{
	    
	    if($tok_len > DEFAULT_LINE_LENGTH)
	    {
		if($line_beg)
		{
		    $new_string .= $tok . "\n" . $char_to_start . $indent;
		    $line_beg = 1;
		    $curr_line_len = $indent_len;
		}
		else
		{
		    $new_string .= "\n" . $char_to_start . $indent . $tok . "\n" . $char_to_start . $indent;
		    $line_beg = 1;
		    $curr_line_len = $indent_len;
		}
	    }
	    elsif(($curr_line_len + $tok_len) > DEFAULT_LINE_LENGTH)
	    {
		#print("adding \\n\n");
		#print("tok is $tok\n");
		#print("curr line len is $curr_line_len\n");
		#print("tok len is $tok_len\n");
		#print("indent_len is $indent_len\n");
		$new_string .= "\n" . $char_to_start . $indent . $tok . " ";
		$line_beg = 1;
		$curr_line_len = $indent_len + $tok_len + 1;
	    }
	    else
	    {
		#print("tok is $tok\n");
		$new_string .= $tok . " ";
		$curr_line_len += ($tok_len + 1);
		#print("curr_line_len now $curr_line_len\n");
		$line_beg = 0;
	    }
	}
    }
    return $new_string;
}

#################################################################
# subroutine : create_special_line
# sub class  : general 
#
# EPN 03.03.05
# 
# purpose : Given a character, return a line that is simply
#           that character repeated $length times
#           followed by a "\n"
# 
# args : (1) $char
#            special character
#        (2) $length
#            number of times to repeat (1)
#################################################################
sub create_special_line
{
    my ($char, $length) = @_;
    my $i;

    my $spec_line = "";
    for($i = 0; $i < $length; $i++)
    {
	$spec_line .= $char;
    }
    $spec_line .= "\n";
    return $spec_line;
}

#################################################################
# subroutine : array_to_scalar
# sub class  : general
# 
# EPN 07.04.05 Driving from Prophetstown, IL --> St. Louis
# 
# purpose : Concatentate all the members of an array into a
#           scalar and return that value.
#
# args : (1) $arr_ref 
#            reference to the array to concatenate
################################################################# 

sub array_to_scalar
{
    my ($arr_ref) = $_[0];
    my $to_return = "";
    my $el;
    foreach $el (@{$arr_ref})
    {
	$to_return .= $el;
    }
    return $to_return;
}
