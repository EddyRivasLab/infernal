#!/usr/bin/perl
#
# Eric Nawrocki
# 09.12.05
# 
# hmm_generate_embed.pl - This script reads in a parameter file which describes an HMM.
#                         It generates a specified number of 'chromosomes' 
#                         using the HMM, and embeds sequences from a set of different
#                         files randomly throughout the generated sequence.
#           
# Modified from Algorithms for Computational Biology Spring 04 HW4 part 1
# Author: Eric Nawrocki
#############################################################################
use strict;
use constant FASTA_LINE_LENGTH       => 50;
use constant SPECIAL_LINE_LENGTH     => 75; 
use constant DEFAULT_LINE_LENGTH     => 75;   

my $usage ="perl hmm_generate_embed.pl\n\t<params file>\n\t<num output seqs>\n\t<length of output seqs>\n\t<output root>\n\t<index file with names of roots>\n\t<directory with each root.test files with seqs to embed>\n\t<seed for RNG>\n";

(@ARGV == 7) || die $usage;

my($params_file, $num_out_seqs, $output_len, $out_root, $idx_file, $embed_dir, $seed) = @ARGV;

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

my %embed_seq_hash_of_hash; #2D hash of seqs to embed key 1:fam name key 2:seq name
my %all_embed_seq_hash; #seq hash of sequences to embed 
my %fam_by_seq_hash; #key is sequence name; value is family name
my $fam;
my $root;
my $key;
my $embed_file;

#read in the parameters from the file
read_params($params_file, \%init_p_hash, \%trans_p_hash, \%emit_p_hash, \@st_names_arr, \@out_alph_arr);

#produce output based on HMM
#produce_output($output_len, \%init_p_hash, \%trans_p_hash, \%emit_p_hash, \@st_names_arr, \@out_alph_arr, \@output_arr, $num_out_seqs, $out_root);

#for each root in the $idx file, read in the embedded sequences
open (INDEX,$idx_file) || die;
while (<INDEX>) 
{
    if(/^(\S+)/) 
    {
	$root = $1;
	#read the embedded sequences in
	$embed_file = "$embed_dir/$root.test";
	read_fasta($embed_file, \%{$embed_seq_hash_of_hash{$root}});
    }
}
#combine all the embedded sequence hashes into a single hash
%all_embed_seq_hash = ();
foreach $fam (keys(%embed_seq_hash_of_hash))
{
    foreach $key (keys(%{$embed_seq_hash_of_hash{$fam}}))
    {
	$all_embed_seq_hash{$key} = $embed_seq_hash_of_hash{$fam}{$key};
	$fam_by_seq_hash{$key} = $fam;
    }
}
my $out_fa_file = $out_root . ".fa";

#produce output based on HMM with sequences embedded
generate_and_embed($output_len, \%init_p_hash, \%trans_p_hash, \%emit_p_hash, \@st_names_arr, \@out_alph_arr, \@output_arr, $num_out_seqs, $out_root, \%all_embed_seq_hash, \%fam_by_seq_hash, $out_fa_file);

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
    for($i=0;$i<$num_states;$i++)
    {
	$init_p_hash_ref->{$st_names_arr_ref->[$i]} = $init_p_arr[$i];
    }
    #read in trans probs
    for($i=0;$i<$num_states;$i++)
    {
	@trans_p_arr = split " ", <PFILE>;
	for($j=0;$j<$num_states;$j++)
	{
	    $trans_p_hash_ref->{$st_names_arr_ref->[$i]}{$st_names_arr_ref->[$j]} = $trans_p_arr[$j];
	}
    }
    @{$out_alph_arr_ref} = split " ", <PFILE>;
    my $len_alph = scalar(@{$out_alph_arr_ref});
    for($i=0;$i<$num_states;$i++)
    {
	@emit_p_arr = split " ", <PFILE>;
	for($j=0;$j<$len_alph;$j++)
	{
	    $emit_p_hash_ref->{$st_names_arr_ref->[$i]}{$out_alph_arr_ref->[$j]} = $emit_p_arr[$j];
	}
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

sub generate_and_embed
{
    #Given parameters of an HMM and sequences to embed,
    #it randomly produces background sequences from the HMM
    #and randomly chooses spots to embed each of the sequences to 
    #embed.  Ensures no two embedded sequences overlap in a hacky
    #manner.
    #
    #calls choose_from_dist iteratively
    my($output_len, $init_p_hash_ref, $trans_p_hash_ref, $emit_p_hash_ref, 
       $st_names_arr_ref, $out_alph_arr_ref, $output_arr_ref, $num_seqs,
       $name_root, $embed_seq_hash_ref, $embed_fam_by_seq_hash_ref, $out_fa_file) = @_;
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
    my $num_embed_seqs = scalar(keys(%{$embed_seq_hash_ref}));
    my $curr_len;
    my $curr_chrom;
    my $curr_begin;
    my $curr_end;
    my $other_begin;
    my $other_end;
    my $curr_orient;
    my $conflict;
    my $key;
    my @sorted_embed_keys = sort (keys(%{$embed_seq_hash_ref}));
    my %curr_embed_info_hash = ();
    my %chrom_begin_hash = ();
    my %chrom_end_hash = ();
    my $chrom_num;
    my $curr_chrom_num_embed;

    my $curr_state;

    my %by_chrom_hash_of_arr_of_hash = ();
    # D1 : key 1 : chromosome number
    # D2 : array : num elements is number of embedded seqs in this chromosome
    # D3 : 4 element hash : "key"   = embedded seq name
    #                       "begin" = begin posn of embedded seq
    #                       "end"   = end posn of embedded seq
    #                       "orient"= orientation of embedded seq (0 = forward; 1 = rev comp)

    my @curr_sorted_chrom_arr_of_hash = ();
    my %hash = ();
    my @temp_sorted_keys = ();

    my $embed_left;
    my $embed_idx;
    my $next_embed_key;
    my $next_embed_begin;
    my $next_embed_end;
    my $next_embed_orient;
    my $offset;
    my $seq_to_embed;
    my @seq_arr_to_embed;
    my %to_print_seq_hash = ();
    my @ordered_names = ();
    my $embed_out = $out_root . ".ebd";
    my $chr_embed_out;

    my %chrom_to_print_seq_hash = ();
    my @chrom_names = ();
    my $chrom_out_fa_file;
    my $conflict_ctr;

    srand(($seed * 100000)); #time|$$ for "random" seed

    #determine the positions to embed the sequences
    
    #initialize the chromosome begin and end hashes
    for($a=0; $a < $num_seqs; $a++)
    {
	@{$chrom_begin_hash{$a}} = ();
	@{$chrom_end_hash{$a}} = ();
	@{$by_chrom_hash_of_arr_of_hash{$a}} = ();
    }

    foreach $key (@sorted_embed_keys)
    {
	$conflict = 1;
	$curr_len = length($embed_seq_hash_ref->{$key});
	$curr_chrom = int(rand(scalar($num_seqs)));
	$curr_orient = int(rand(2));
	#printf("embedding $key | curr chrom : $curr_chrom | curr_orient : $curr_orient\n");
	if($curr_len > $output_len)
	{
	    printf("ERROR: trying to embed seq $key of length $curr_len into chromosome of length $output_len.\n");
	    exit();
	}
	#printf("\n");
	$conflict_ctr = 0;
	while($conflict)
	{
	    $conflict_ctr++;
	    if($conflict_ctr > 10000)
	    {
		print("ERROR, having trouble embedding a sequence such that it doesn't overlap another embedded sequence.\nTry increasing the size of the sequences or the length of the sequences.\n");
		die;
	    }
	    #we could pick a new chromosome and orientation after each conflict:
	    #$curr_chrom = int(rand(scalar($num_seqs)));
	    #$curr_orient = int(rand(2));
	    $curr_begin = int(rand(($output_len - $curr_len)));
	    $curr_end = $curr_begin + $curr_len - 1;
	    #printf("key : $key | curr_chrom : $curr_chrom | curr_orient : $curr_orient | curr begin : $curr_begin | curr_end : $curr_end\n");	    
	    #check to make sure there are no conflicts with other embedded
	    #sequences in this chromosome.
	    
	    $conflict = 0;
	    for($i = 0; $i < scalar(@{$chrom_begin_hash{$curr_chrom}}); $i++)
	    {
		$other_begin = $chrom_begin_hash{$curr_chrom}[$i];
		$other_end = $chrom_end_hash{$curr_chrom}[$i];
		#printf("checking against $other_begin -> $other_end\n");
		# check for this situation:
		#            (other depicted by .'s)
		#            ..................
		#            __________________
		#            (curr_begin= somewhere within _'s)
		if(($curr_begin >= $other_begin) && ($curr_begin <= $other_end))
		{
		    #print("found conflict with $other_begin -> $other_end\n");
		    $conflict = 1;
		    last;
		}
		# check for this situation:
		#            (other depicted by .'s)
		#            ..................
		#            ++++++++++++++++++
		#            (curr_end= somewhere within +'s)
		if(($curr_end >= $other_begin) && ($curr_end <= $other_end))
		{
		    #print("found conflict with $other_begin -> $other_end\n");
		    $conflict = 1;
		    last;
		}
		# check for this situation:
		#            (other depicted by .'s)
		#            ..................
		#   <--______                  ++++++-->  
		#            (curr_begin = somewhere within _'s (or to left of _'s)
		#            (curr_end   = somewhere within +'s (or to right of +'s)
		if(($curr_begin < $other_begin) && ($curr_end > $other_end))
		{
		    #print("found conflict with $other_begin -> $other_end\n");
		    $conflict = 1;
		    last;
		}
	    }
	    #now if conflict = 0 ==> there are no conflicts
	}
	
	%curr_embed_info_hash = ();
	$curr_embed_info_hash{"key"} = $key;
	$curr_embed_info_hash{"begin"} = $curr_begin;
	$curr_embed_info_hash{"end"} = $curr_end;
	$curr_embed_info_hash{"orient"} = $curr_orient;

	$curr_chrom_num_embed = scalar(@{$by_chrom_hash_of_arr_of_hash{$curr_chrom}});
	%{$by_chrom_hash_of_arr_of_hash{$curr_chrom}[$curr_chrom_num_embed]} = %curr_embed_info_hash;
	push(@{$chrom_begin_hash{$curr_chrom}}, $curr_begin);
	push(@{$chrom_end_hash{$curr_chrom}}, $curr_end);
    }

    #print out information about each embedded sequence 
    open(OUT, ">" . $embed_out);
    foreach $chrom_num (sort keys (%by_chrom_hash_of_arr_of_hash))
    {
	$chr_embed_out = $embed_out;
	$chr_embed_out =~ s/\.ebd$//;
	$chr_embed_out .= "_chr" . ($chrom_num+1) . ".ebd";

	open(CHROUT, ">" . $chr_embed_out);
	for($i = 0; $i < (scalar(@{$by_chrom_hash_of_arr_of_hash{$chrom_num}})); $i++)
	{
	    printf OUT ("%s %s %s %d %d %d\n", 
			$embed_fam_by_seq_hash_ref->{$by_chrom_hash_of_arr_of_hash{$chrom_num}[$i]{"key"}},
			$by_chrom_hash_of_arr_of_hash{$chrom_num}[$i]{"key"},
			($name_root . "_" . ($chrom_num+1)),
			($by_chrom_hash_of_arr_of_hash{$chrom_num}[$i]{"begin"}+1),
			($by_chrom_hash_of_arr_of_hash{$chrom_num}[$i]{"end"}+1),
			$by_chrom_hash_of_arr_of_hash{$chrom_num}[$i]{"orient"});
	    printf CHROUT ("%s %s %s %d %d %d\n", 
			   $embed_fam_by_seq_hash_ref->{$by_chrom_hash_of_arr_of_hash{$chrom_num}[$i]{"key"}},
			   $by_chrom_hash_of_arr_of_hash{$chrom_num}[$i]{"key"},
			   ($name_root . "_" . ($chrom_num+1)),
			   ($by_chrom_hash_of_arr_of_hash{$chrom_num}[$i]{"begin"}+1),
			   ($by_chrom_hash_of_arr_of_hash{$chrom_num}[$i]{"end"}+1),
			   $by_chrom_hash_of_arr_of_hash{$chrom_num}[$i]{"orient"});
	    #debug_print_hash(\%{$by_chrom_hash_of_arr_of_hash{$chrom_num}[$i]}, "embed seq num $i\n");
	}
	close(CHROUT);
    }
    close(OUT);

    print_out_file_notice($embed_out, "File containing information on $num_embed_seqs embedded sequences and where they appear in the file $out_fa_file. Format per line: <s1> <s2> <s3> <d1> <d2> <d3>. <s1> = name of family embedded sequence belongs to. <s2> = name of embedded sequence. <s3> = name of sequence its embedded in. <d1> = beginning position of embedded seq (index starting from 1) <d2> = end position of embedded sequence. <d3> = 0 or 1 for either forward or reverse orientation.\n Each chromosomes embedded sequences were printed to *chr*.ebd files also.\n");
    
    #produce output of specified length based on the HMM parameters
    $curr_state = $st_names_arr_ref->[rand(scalar (keys %{$init_p_hash_ref}))];

    #for each sequence (chromosome in case of pseudo-genome)
    for($i=0; $i < $num_seqs; $i++)
    {
	$name = $name_root . "_" . ($i + 1);
	#printf("\n\nchrom: " . $name . "\n");

	#create the output 
	@{$output_arr_ref} = ();
	
	#first we have to sort our embedded sequences by start positions
	#for the current output sequence (chromosome)
	@curr_sorted_chrom_arr_of_hash = ();
	%hash = ();
	for($k = 0; $k < scalar(@{$by_chrom_hash_of_arr_of_hash{$i}}); $k++)
	{
	    $hash{$k} = $by_chrom_hash_of_arr_of_hash{$i}[$k]{"begin"};
	}
	#@temp_sorted_keys = sort sort_hash_value_num_inc keys (%hash);
	@temp_sorted_keys = sort { $hash{$a} <=> $hash{$b} } keys (%hash);
	#print("\n");
	#debug_print_arr(\@temp_sorted_keys, "temp sorted keys");
	#print("\n");
	for($k = 0; $k < scalar(@temp_sorted_keys); $k++)
	{
	    %{$curr_sorted_chrom_arr_of_hash[$k]} = %{$by_chrom_hash_of_arr_of_hash{$i}[$temp_sorted_keys[$k]]};
	}
	
	if(scalar(@curr_sorted_chrom_arr_of_hash) == 0)
	{
	    $embed_left = 0;
	    $embed_idx = -1;
	    $next_embed_key = -1;
	    $next_embed_begin = -1;
	    $next_embed_end = -1;
	    $next_embed_orient = -1;
	}
	else
	{
	    $embed_left = 1;
	    $embed_idx = 0;
	    $next_embed_key = $curr_sorted_chrom_arr_of_hash[0]{"key"};
	    $next_embed_begin = $curr_sorted_chrom_arr_of_hash[0]{"begin"};
	    $next_embed_end = $curr_sorted_chrom_arr_of_hash[0]{"end"};
	    $next_embed_orient = $curr_sorted_chrom_arr_of_hash[0]{"orient"};
	}
	for($q=0;$q<$output_len;$q++)
	{
	    #check if we need to embed, if so, we skip i ahead
	    if($q == $next_embed_begin)
	    {
		#printf("q : $q embedding seq : $next_embed_key\n");
		if($next_embed_orient)
		{
		    $seq_to_embed = rev_comp_rna($embed_seq_hash_ref->{$next_embed_key});
		}
		else
		{
		    $seq_to_embed = $embed_seq_hash_ref->{$next_embed_key};
		}
		@seq_arr_to_embed = split("", $seq_to_embed);
		$offset = $q;
		for(;$q <= $next_embed_end; $q++)
		{
		    $output_arr_ref->[$q] = $seq_arr_to_embed[($q-$offset)];
		}
		# determine attributes of next seq to embed (if there's any left)
		if(scalar(@curr_sorted_chrom_arr_of_hash) == ($embed_idx+1))
		{
		    #print("embed not left\n");
		    $embed_left = 0;
		    $embed_idx = -1;
		    $next_embed_key = -1;
		    $next_embed_begin = -1;
		    $next_embed_end = -1;
		    $next_embed_orient = -1;
		}
		else
		{
		    #print("embed left\n");
		    $embed_left = 1;
		    $embed_idx++;
		    $next_embed_key = $curr_sorted_chrom_arr_of_hash[$embed_idx]{"key"};
		    $next_embed_begin = $curr_sorted_chrom_arr_of_hash[$embed_idx]{"begin"};
		    $next_embed_end = $curr_sorted_chrom_arr_of_hash[$embed_idx]{"end"};
		    $next_embed_orient = $curr_sorted_chrom_arr_of_hash[$embed_idx]{"orient"};
		    #print("next embed begin : $next_embed_begin\n");
		}
	    }

	    #when we get here, the HMM is in charge of emitting chromosome posn $q
	    #its possible that our embedded seq took up the rest of the output so
	    #we have to check
	    #or our next embedded seq might be immediately after our last one ended
	    #we have to check for this too.
	    if(($q < $output_len) && ((!($embed_left)) || ($q < $next_embed_begin)))
	    {
		#
		#step 1 emit letter of output alphabet based on current state
		#create probability dist for emission of current state
		$j=0;
		foreach $char (@{$out_alph_arr_ref})
		{
		    $curr_prob_dist_arr[$j++] = $emit_p_hash_ref->{$curr_state}{$char};
		}
		#emit letter of output alphabet based on probability dist
		$output_arr_ref->[$q] = $out_alph_arr_ref->[choose_from_dist(\@curr_prob_dist_arr)];
		
		#step 2, pick a new state based on transition probabilities
		#create probability dist for transition from current state
		@curr_prob_dist_arr = ();
		$j=0;
		foreach $state (@{$st_names_arr_ref})
		{
		    @curr_prob_dist_arr[$j++] = $trans_p_hash_ref->{$curr_state}{$state};
		}
		#transition to new state based on probability dist
		$curr_state = $st_names_arr_ref->[choose_from_dist(\@curr_prob_dist_arr)];
	    }
	    #hack to ensure that if two embedded seqs are right next to each other
	    #then we have to decrement i so the for loop doesn't increment i past
	    #the beginning of the next embedded sequence
	    if(($embed_left) && ($q == $next_embed_begin))
	    {
		$q--;
	    }
	}
	#print output
	#for($i=0;$i<$output_len;$i++)
	#{
	#    if(($i != 0) && ($i % 60 == 0)) { print("\n"); }
	#    #print $output_arr_ref->[$i];
	#}
	#print("pushed name to ordered names\n");

	#print out this chromosome to its own fasta file
	%chrom_to_print_seq_hash = ();
	$chrom_to_print_seq_hash{$name} = array_to_scalar($output_arr_ref);
	@chrom_names = ($name);
	$chrom_out_fa_file = $out_fa_file;
	$chrom_out_fa_file =~ s/\.fa//;
	$chrom_out_fa_file .= "_chr" . ($i+1) . ".fa";
	print_fasta_ordered(\%chrom_to_print_seq_hash, $chrom_out_fa_file, \@chrom_names);

	$to_print_seq_hash{$name} = array_to_scalar($output_arr_ref);
	push(@ordered_names, $name);
    }
    #print out the whole genome to a fasta file
    print_fasta_ordered(\%to_print_seq_hash, $out_fa_file, \@ordered_names);
    print_out_file_notice($out_fa_file, "New fasta file containing the " . (scalar(@ordered_names)) . " sequences generated by HMM explained in $params_file, and with the $num_embed_seqs embedded sequences randomly inserted from file $embed_file.");
    
    my $chrom_files = $out_fa_file;
    $chrom_files =~ s/\.fa//;
    my $chrom_files_begin = $chrom_files . "_chr1.fa";
    my $chrom_files_end = $chrom_files . "_chr" . scalar(@ordered_names) . ".fa";
    $chrom_files = $chrom_files_begin . "->" . $chrom_files_end;

    print_out_file_notice($chrom_files, "Each chromosome in $out_fa_file in its own fasta file.");
    
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


#################################################################
# subroutine : read_fasta_careful
# sub class  : crw and sequence
# 
# EPN 01.26.06
#
# purpose : Open, read, and store the information in a given
#           .fa (fasta format) file. CAREFUL because it 
#           removes everything in the header starting at
#           the first space.
#
# args : (1) $in_file
#            name of .fa file in current directory
#        (2) $seq_hash_ref
#            reference to the hash that will contain the sequence
#            information.  Fasta description line used as key for
#            each sequence, sequence is value.
################################################################# 
sub read_fasta_careful
{
    my ($in_file, $seq_hash_ref) = @_;
    open(IN, $in_file);

    my $line;
    my $seq_name;

    #chomp up beginning blank lines
    $line = <IN>;
    while(!($line =~ m/^>/))
    {
	$line = <IN>;
    }

    chomp $line;
    $line =~ s/\s.*//;
    $seq_name = $line;
    $seq_name =~ s/^>//;
    while($line = <IN>)
    {
	chomp $line;
	$seq_hash_ref->{$seq_name} = "";
	while((!($line =~ m/^>/)) && ($line ne ""))
	{
	    $seq_hash_ref->{$seq_name} .= $line;
	    $line = <IN>;
	    chomp $line;
	}
	chomp $line;
	$line =~ s/\s.*//;
	$seq_name = $line;
	$seq_name =~ s/^>//;
    }
    #trim_keys_in_hash($seq_hash_ref, DEFAULT_MAX_SEQ_HEADER_LENGTH);
}
