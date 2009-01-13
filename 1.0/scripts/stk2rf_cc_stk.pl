#!/usr/bin/perl
#
#
# Eric Nawrocki
# EPN, Wed Feb  7 09:54:12 2007
#
# stk2rf_cc_stk.pl Given a stockholm (.stk) alignment file and a gap threshold.
#                  add a #=GR RF line with x's and .'s with a x at each position
#                  with less than the gap threshold gaps. Also mark the consensus
#                  position of each column with #=GR COL lines (as many as nec
#                  to annotate all columns, for example if there's 100 columns
#                  we need 3 #=GR lines for each column.)
#
# Subroutines copied from M_seq.pm and M_gen.pm.
#                   
#
# General Strategy
#
# (A) Read in original stk alignment file
# (B) Determine number of gaps in each column
# (C) Create and print new stockholm file
use constant DEFAULT_LINE_LENGTH     => 75;   
use constant SPECIAL_LINE_LENGTH     => 75; 
use constant STK_LINE_LENGTH         => 50;

$usage = "perl stk2rf_stk.pl <stk file> <gap threshold> <output stk file name>\n";

if(@ARGV != 3)
{
    die $usage;
}

($in_file, $gap_thresh, $out_file) = @ARGV;

#####################################################################
# (A) Read in original stk alignment file
#####################################################################
%full_seq_hash = ();
%bp_hash = ();

#get a seq_hash and bp_hash (which is unused) from the alignment
read_stk_all_markup_as_seqs($in_file, \%full_seq_hash, 1);
%seq_hash = ();
foreach $key (keys(%full_seq_hash))
{
    if($key !~ m/^\#/)
    {
	$seq_hash{$key} = $full_seq_hash{$key};
    }
}

#first get the columns into a manageable data structure
@cols_arr;
get_cols_from_aln(\%seq_hash, \@cols_arr);

#determine which columns have gaps above the threshold
@gap_status;
characterize_aln_cols_gap_thresh(\@cols_arr, \@gap_status, $gap_thresh);

$rf_string = "";
$cc = 1;
for($i = 1; $i < scalar(@gap_status); $i++)
{
    #printf("gap_status: $i $gap_status[$i]\n");
    if($gap_status[$i] == 0)
    {
	$rf_string .= ".";
	$cons_col_A[$i] = 0;
    }
    else
    {
	$rf_string .= "x";
	$cons_col_A[$i] = $cc++;
    }
}
#printf("rf string: $rf_string\n");
$full_seq_hash{"#=GC RF"} = $rf_string;
#Now create the '#=GC COL' markup 
$tmp_cc = $cc;
$ndigits = 1;
while($tmp_cc > 10)
{
    $tmp_cc = int($tmp_cc / 10);
    $ndigits++;
}
for($i = 1; $i <= $ndigits; $i++)
{
    $new_key = "#=GC COL ";
    for($j = 1; $j < $i; $j++)
    {
	$new_key .= "x";
    }
    $new_key .= "A";
    for($j = $ndigits; $j > $i; $j--)
    {
	$new_key .= "x";
    }
    $new_markup = "";

    $nzeroes = $ndigits - $i;
    $denominator = 1;
    for($x = 1; $x <= $nzeroes; $x++)
    {
	$denominator *= 10;
    }
    for($x = 1; $x < scalar(@gap_status); $x++)
    {
	if($cons_col_A[$x] != 0)
	{
	    $curr_val = int($cons_col_A[$x] / $denominator);
	    $curr_char= substr($curr_val, -1, 1);
	    #printf("cc: $cons_col_A[$x] denom: $denominator curr_val: $curr_val, curr_char: $curr_char\n");
	    
	    $new_markup .= $curr_char
	}
	else
	{
	    $new_markup .= ".";
	}
    }	
    $full_seq_hash{$new_key} = $new_markup;
}

print_stk_no_struct(\%full_seq_hash, $out_file, "Same aln as $in_file, RF line deduced with gap threshold of $gap_thresh. COL lines denote consensus column positions.");

print_out_file_notice($out_file, "Same alignment as $in_file with a RF line added deduced using a gap threshold of $gap_thresh. COL lines denote consensus column postions.");

#################################################################################
#SUBROUTINES: copied from perl modules
#################################################################################


#################################################################
# subroutine : read_stk_all_markup_as_seqs
# sub class  : sequence 
# 
# EPN 06.13.05
#
# purpose : Open, read, and store the information in a given
#           .stk (stockholm format) file.  Treat all markup
#           #=GR and #=GC lines as individual sequences, that is
#           give them each
#           their own value in the seq hash. (For example : 
#           #=GC SS_cons lines would be stored in the hash
#           key "#=GC SS_cons" where the value would be the
#           secondary structure sequence.
#
# args : (1) $in_file
#            name of .stk file in current directory
#        (2) $seq_hash_ref
#            reference to the hash that will contain the sequence
#            information.  Name of each seq will be used as key for
#            each sequence, sequence is value.
################################################################# 
sub read_stk_all_markup_as_seqs
{
    ($in_file, $seq_hash_ref) = @_;
    open(IN, $in_file);
    
    #print("in file $in_file\n");

    $line = <IN>;
    #we don't care about any commented lines except #=
    # - this is probably not a good
    #idea at all

    #assumes first line we care about 
    #is a regular (non-markup) sequence line
    while(($line =~ m/^\#/) || ($line !~ m/\w/))
    {
	$line = <IN>;
	#print("read line $line\n");
    }

    #get the first seq
    ($curr_name, $curr_seq) = split(/\s+/, $line);
    $seq_hash_ref->{$curr_name} .= $curr_seq;

    while($line = <IN>)
    {
	if($line =~ m/^\w/)
	{
	    #regular (non-markup) sequence line
	    ($curr_name, $curr_seq) = split(/\s+/, $line);
	    $seq_hash_ref->{$curr_name} .= $curr_seq;
	}
	elsif($line =~ m/^\#=GR/)
	{
	    #GR markup line
	    ($GR, $curr_name, $curr_markup, $curr_seq) = split(/\s+/, $line);
	    $curr_key = $GR . " " . $curr_name . " " . $curr_markup;
	    $seq_hash_ref->{$curr_key} .= $curr_seq;
	}
	elsif($line =~ m/^\#=GC/)
	{
	    #GC markup line
	    ($GC, $curr_name, $curr_seq) = split(/\s+/, $line);
	    $curr_key = $GC . " " . $curr_name;
	    $seq_hash_ref->{$curr_key} .= $curr_seq;
	}
    }
    #trim_keys_in_hash($seq_hash_ref, DEFAULT_MAX_SEQ_HEADER_LENGTH);
    #trim_keys_in_hash($mark_up_seq_hash_of_hash_ref, DEFAULT_MAX_SEQ_HEADER_LENGTH);
}

################################################################
# subroutine : get_cols_from_aln
# sub class  : sequence
# 
# EPN 03.08.05
#
# purpose : Given an alignment fill an array with each column
#           of the alignment as an array for the value of that
#           hash and the key is the column of the alignment.
#
# args : (1) $seq_aln_hash_ref
#            reference to a hash that is the alignment
#        (2) $cols_arr_ref
#            reference to an array, for which posn[x] will
#            be filled with a string containgin each residue in
#            posn x of the alignment.
#################################################################
sub get_cols_from_aln
{
    ($seq_aln_hash_ref, $cols_arr_ref) = @_;
    $aln_length = get_aln_length($seq_aln_hash_ref);
    
    $cols_arr_ref->[0] = "POSN 0 IS BLANK, SEQ INDEXING STARTS AT 1";

    for($j = 0; $j < $aln_length; $j++)
    {
	$position = $j + 1;
	foreach $header (sort keys(%{$seq_aln_hash_ref}))
	{
	    $cols_arr_ref->[$position] .= 
		substr($seq_aln_hash_ref->{$header}, $j, 1)
	}
    }
}

#################################################################
# subroutine : characterize_aln_cols_gap_thresh
# sub class  : sequence
#
# EPN 04.30.06
# (derived from characterize_ae2_gb_aln_cols sub in M_crw.pm)
# 
# purpose : Given an array of strings that are each a column of 
#           an alignment, determine which columns have a fraction
#           of sequences that are gaps above a given threshold.
#           (all non-word characters), and save that info
#           by setting col_status_arr_ref->[i] to 0
#           for column i if gap threshold is EXCEEDED
#
#           Treats all non-word characters as gaps, so be 
#           careful to remove all structure annotations etc.
#           prior to making @cols_arr.
# args : (1) $cols_arr_ref
#            reference to an array that is the columns of the 
#            alignment
#        (2) $col_status_arr_ref
#            reference to an array to be filled with 1s or 0s.
#            Position[i] = 1 if we want to keep that column, 0
#            otherwise.
#        (3) $gap_thresh
#            if the fraction of seqs in col i with gaps exceeds this number
#            set col_status_arr_ref->[i] to 0.
#################################################################
sub characterize_aln_cols_gap_thresh
{
    my ($cols_arr_ref, $col_status_arr_ref, $gap_thresh) = @_;

    #we index sequences starting from position 1
    $col_status_arr_ref->[0] = -1;

    for($col_posn = 1; $col_posn < scalar(@{$cols_arr_ref}); $col_posn++)
    {
	$has_res = 0;
	#print("checking column $col_posn\n");
	$aln_col = $cols_arr_ref->[$col_posn];
	$orig_len = length($aln_col);
	#remove gaps (all non-word characters)
	$aln_col =~ s/\W//g;
	$nongap_len = length($aln_col);
	$gap_fract = 1.0 - ($nongap_len / $orig_len);
	#printf("i: $col_posn | gap_fract: $gap_fract| orig_len: $orig_len | nongap_len: $nongap_len\n");
	if($gap_fract > $gap_thresh) 
	{
	    $col_status_arr_ref->[$col_posn] = 0;
	}
	else
	{
	    $col_status_arr_ref->[$col_posn] = 1;
	}	    
    }
}


#################################################################
# subroutine : print_stk_no_struct
# sub class  : sequence
#
# EPN 03.03.05
# 
# purpose : Print a stockholm formatted file with no structure
# 
# args : (1) $seq_hash_ref
#            reference to a sequence hash
#        (2) $out_file_name
#            name of file to print to
#        (3) $extra_commented_line
#            extra line to be added after the #STOCKHOLM 1.0 
#################################################################
sub print_stk_no_struct
{
    ($seq_hash_ref, $out_file_name, $extra_commented_line) = @_;
    
    #determine longest header
    @headers = keys (%{$seq_hash_ref});
    $max_head_len = get_max_len_string_from_arr(\@headers);
    $header_len = $max_head_len + 10;

    $max_seq_len = get_max_len_value_string_from_hash($seq_hash_ref);
    
    open(OUT, ">" . $out_file_name);
    print OUT "# STOCKHOLM 1.0\n";
#    print OUT "#=GF ?\n";
    if($extra_commented_line ne "")
    {
	$new_line = add_newlines_to_string($extra_commented_line, 0, 0, "#=GF AU "); 
	print OUT $new_line . "\n" . "\n";
    }
    else
    {
	print OUT "\n";
    }
    $index = 0;
    while($index < $max_seq_len)
    {
	@last_lines = ();
	foreach $key (sort keys (%{$seq_hash_ref}))
	{
	    #only print it if the key is not a #=GC line
	    if($key =~ m/^\#=GC/)
	    {
		push(@last_lines, $key);
	    }
	    else
	    {
		$curr_substr = substr($seq_hash_ref->{$key}, $index, STK_LINE_LENGTH);
		printf OUT ("%-" . $header_len . "s $curr_substr\n", $key);
	    }
	}
	foreach $key (@last_lines)
	{
	    $curr_substr = substr($seq_hash_ref->{$key}, $index, STK_LINE_LENGTH);
	    printf OUT ("%-" . $header_len . "s $curr_substr\n", $key);
	}	    
	$index += STK_LINE_LENGTH;
	print OUT ("\n");
    }
    print OUT ("//\n");
    close(OUT);
}
#################################################################
# subroutine : get_aln_length
# sub class  : sequence
# 
# EPN 03.08.05
#
# purpose : Given an alignment determine the length of each
#           sequence in the alignment.  If they are not all the
#           same length, die and print an error message.  If they
#           are all the same length, return the length.
#
# args : (1) $aln_hash_ref
#            reference to a hash that is the alignment
#################################################################
sub get_aln_length
{
    $aln_hash_ref = $_[0];
    @sorted_keys =  (sort keys %{$aln_hash_ref});
    $aln_length = length($aln_hash_ref->{$sorted_keys[0]});

    #check to make sure all sequences are same length
    for($i = 1; $i < scalar(@sorted_keys); $i++)
    {
	if($aln_length != length($aln_hash_ref->{$sorted_keys[$i]}))
	{
	    print("ERROR in get_aln_length\n\tnot all sequences in alignment are of the same length\n");
	    print("aln length is $aln_length\n");
	    print("but $sorted_keys[$i] has length " . length($aln_hash_ref->{$sorted_keys[$i]}) . "\n");
	    print("seq is $aln_hash_ref->{$sorted_keys[$i]}\n");
	    die;
	}
    }
    return $aln_length;
}

#################################################################
# subroutine : get_max_len_string_from_arr
# sub class  : general
#
# EPN 03.03.05
# 
# purpose : Given an array of strings, return the length of the
#           longest string in the array
# 
# args : (1) $arr_ref
#            reference to array of strings
#################################################################
sub get_max_len_string_from_arr
{
    $arr_ref = $_[0];

    $max_len = length($arr_ref->[0]);
    for($i = 1; $i < scalar(@{$arr_ref}); $i++)
    {
	if(length($arr_ref->[$i]) > $max_len)
	{
	    $max_len = length($arr_ref->[$i]);
	}
    }
    return $max_len;
}

#################################################################
# subroutine : get_max_len_value_string_from_hash
# sub class  : general
#
# EPN 03.03.05
# 
# purpose : Given a hash with strings as values, return the length of the
#           longest string of all the values
# 
# args : (1) $hash_ref
#            reference to hash with strings as values
#################################################################
sub get_max_len_value_string_from_hash
{
    $hash_ref = $_[0];

    @sorted_keys = sort keys(%{$hash_ref});
    $max_len = length($hash_ref->{$sorted_keys[0]});
    for($i = 1; $i < scalar(@sorted_keys); $i++)
    {
	if(length($hash_ref->{$sorted_keys[$i]}) > $max_len)
	{
	    $max_len = length($hash_ref->{$sorted_keys[$i]});
	}
    }
    return $max_len;
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
    ($string, $indent_len, $first_line_length, $char_to_start) = @_;

    @string_arr = split(" ", $string);
    $line_beg = 1;
    $new_string = $char_to_start;
    $curr_line_len = $first_line_length;
    $indent = "";
    for($i = 0; $i < $indent_len; $i++)
    {
	$indent .= " ";
    }

    foreach $tok (@string_arr)
    {
	$tok_len = length($tok);
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
    return $new_string;
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
    ($file_name, $description) = @_;

    $char = "*";
    $spec_line = create_special_line($char, SPECIAL_LINE_LENGTH);
    print("$spec_line");
    print(" Output file notice\n");
    print(" File name   : $file_name\n");
    
    $new_desc = add_newlines_to_string($description, length(" description : "),
				       length(" description : "));
    print(" description : $new_desc\n");
    print("$spec_line");
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
    ($char, $length) = @_;
    
    $spec_line = "";
    for($i = 0; $i < $length; $i++)
    {
	$spec_line .= $char;
    }
    $spec_line .= "\n";
    return $spec_line;
}
