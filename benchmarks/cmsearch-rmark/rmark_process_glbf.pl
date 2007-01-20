#! /usr/bin/perl
#
# Eric Nawrocki 10.11.05
# rmark_process_glbf.pl
#
# Usage:    perl rmark_process_glbf.pl 
#                <.rmm file used to get *.glbf> 
#                <.rmk used to get *.glbf>
#                <seq directory>
#                <index file with fam names; provide path>
#                <genome root <X>, <X.fa> and <X>.ebd  must be seq dir>
#                <concatenated *.glbf files from >= 1 rmark.pl runs; in CWD>
#                <output root>
#
# Options: 
#        Hit resolution options:
#        -R hit : [default] each hit is a single positive/negative
#        -R fnt : treat every nucleotide as a separate positive or negative.
#        -R nnt : treat every non-positive nucleotide as separate negative,
#                and every positive nt as a 1/length(hit) fraction of a hit
#
#        Ignore cross-hits (hits to fam Y while searching with fam X) options:
#        -I both: [default] ignore cross hits on both strands
#        -I none: don't ignore cross hits on either strand
#        -I opp : don't ignore cross hits on opposite strand
#
# Example:  perl rmark_process_glbf.pl infernal.rmm inf-72.rmk 
#                rmark-test/ rmark-test/rmark-test.idx rmark-test 
#                rmark-test_out.glbf rmark-test_out
#
# The example run above will result in the following files:
# rmark_test_out.fam: sorted list of positives and negatives per family along with
#                     stats (MER, >noise, etc.) per family from the index file with fam names
# rmark_test_out.all: sorted list of all positives and negatives along with 
#                     stats on all hits (across all families)
# rmark_test_out.roc: list of x,y points for an ROC curve FOR ALL FAMILIES in format:
#                     '(x,y)' per line
#
#########################################################################################
# Explanation of the options:
#########################################################################################
# The first two options change the hit resolution mode:
#        -F: treat every nucleotide as a separate hit.
#        -N: treat every nucleotide as a 1/length(hit) fraction of a hit
#
# There are 3 different hit resolution modes:
# (1) hit: this is default 
# (2) fnt: enabled with the -F option
# (3) nnt: enabled with the -N option
#
# The differences of these 3 modes is based on how they classify positives and negatives:
# 
# Classification of positives:
#
# If run in 'fnt' mode enabled by the : each nt of the true test ncRNAs (positive sequences) is considered
# a separate positive. The number of positives + negatives (see 'Classification of negatives:'
# below) per family is equal to 2 * the genome length (due to 2 orientations) in this mode.
# So for a benchmark with 51 families and 450 ncRNAs spanning 101,855 nt of a 1 MB pseudo-genome, 
# there are 101,855 positives, and 51 * 2 * 1,000,000 - 101,855 = 
# IF, however the ignore cross hits option is set to "B" (the eigth command line argument) 
# this number of negatives 51 * 2 * (1,000,000 - 101,855) = 91,610,790 negatives (see 'IGNORE
# CROSS HITS OPTION' below).
# AND IF, however the ignore cross hits option is set to "Y" 
# this number of negatives 51 * ((2 * 1,000,000) - 101,855) = 96,805,395 negatives (see 'IGNORE
# CROSS HITS OPTION' below).
#
# If run in 'nnt' mode: each nt of the true test ncRNAs (positive sequences) is considered
# a separate positive, but the weight of a hit is determined by the length of the sequence, each
# nucleotide counts as 1/len(pos) of a positive sequence. This normalizes for length of different
# ncRNAs. Without this normalization a method that always detects the longer ncRNAs will seem
# its better than it really is if we think the importance of detecting the ncRNA length-independent.
#
# If run in 'hit' mode: a hit returned by the search procedure is called a positive if
# it overlaps with an embedded true test sequence and the region of overlap between the
# hit and that true test sequence is > 50% the length of the shorter of the hit or the true
# test sequence. Only one hit to each true test sequence, the best hit as determined by the score, 
# is kept. Nucleotides from hits determined to be 'positives' in this way that overhang into
# 'negative space' (go past the boundary of the true test sequence in the genome) are ignored.
# 
# Classification of negatives:
#
# If run in 'fnt' or 'nnt' mode, each nt of the test sequence that is not part of an embedded
# positive is considered to be a separate 'negative'.
# In these modes, for a chromosome of length 10K with 1 positive sequence of length 100nt embedded 
# in it, there are 19,900 true negatives (one for each negative nt in each orientation).
# For each nucleotide, the script keeps track of the best scoring hit that includes 
# that nucleotide (importantly in nnt mode, only negative hits are considered here, while
# in fnt mode, a positive hit that bleeds into negative sequence results in negative hits in the 
# bleed over region that have the score of the positive hit). If a nucleotide is hit by 
# zero hits using the search method, it is implicitly treated as having the WORST possible score
# for that method.
# 
# If run in 'hit' mode: each entire hit returned by the search method is a separate negative.
# No two negative hits that overlap significantly are allowed. If a negative hit does overlap 
# significantly with another, the better scoring one is kept, and the lower scoring one is
# ignored. Two hits 'significantly overlap' if the region of overlap between them is > 50% the 
# length of the shorter of the two hits.
#
#########################################################################################
# The second two options change how the script treat 'cross-hits'
#        -I: don't ignore cross hits (hits to fam Y with fam X)
#        -O: don't ignore cross hits on other strand
#
# Decided to ignore cross hits after the embarassing realization that all 17 
# of the negatives for RF00169 for glocal cmsearch were SRP_euk_arch sequences,
# RF00169 is SRP_bact. Derr...
#
# Default behavior is, when searching with fam X, to ignore hits that significantly
# overlap with EITHER STRAND of true positives from familiy Y (with X!=Y).
# (doesn't print them AT ALL, and ignores them in MER statistic calculations). 
#
# Enabling the -I options doesn't ignore any cross-hits.
# Enabling the -O option only ignores hits to family Y if they're on the opposite
# strand of the true positive in family Y.
#
#########################################################################################

use Getopt::Std;
$res_opt    = "hit";
$ignore_opt = "both";

getopts('R:I:');
if (defined $opt_R) { $res_opt    = $opt_R; }
if (defined $opt_I) { $ignore_opt = $opt_I; }

$usage = "Usage: perl rmark_process_glbf.pl\n\t<.rmm file used>\n\t<.rmk file used>\n\t<seq directory with *.ali, *.test, *.idx, *.raw files>\n\t<index file with family names; provide path>\n\t<genome root <X>, <X>.fa and <X>.ebd must be in seq dir>\n\t<concatenated *.glbf output from >= 1 rmark.pl runs; in CWD>\n\t<output root>\n";
$options_usage  = "\nOptions: (see code for details)\n\t";
$options_usage .= "Hit resolution options:\n\t";
$options_usage .= "-R hit : [default] each hit is a single positive/negative\n\t";
$options_usage .= "-R fnt : treat every nucleotide as a separate positive or negative.\n\t";
$options_usage .= "-R nnt : treat every non-positive nucleotide as separate negative,\n\t";
$options_usage .= "        and every positive nt as a 1/length(hit) fraction of a hit\n\n\t";
$options_usage .= "Ignore cross-hits (hits to fam Y while searching with fam X) options:\n\t";
$options_usage .= "-I both: [default] ignore cross hits on both strands\n\t";
$options_usage .= "-I none: don't ignore cross hits on either strand\n\t";
$options_usage .= "-I opp : don't ignore cross hits on opposite strand\n";

if(@ARGV != 7)
{
    print $usage;
    print $options_usage;
    exit();
}

$rmm = shift;
$rmk = shift;
$dir = shift;
$idx = shift;
$genome_root = shift;
$glbf_file = shift;
$out_root = shift;

$genome_file = $dir . "/" . $genome_root . ".fa";
$embed_file = $dir . "/" . $genome_root . ".ebd";

if($res_opt eq "hit")
{
    $fnt = 0;
    $nnt = 0;
    $hit = 1;
}
elsif($res_opt eq "nnt")
{
    $fnt = 0;
    $nnt = 1;
    $hit = 0;
}
elsif($res_opt eq "fnt")
{
    $fnt = 1;
    $nnt = 0;
    $hit = 0;
}
else
{
    print ("INVALID -R option\n");
    print $usage;
    print $options_usage;
    exit();
}
if($ignore_opt eq "opp")
{
    $ignore_flag = "opposite";
}
elsif($ignore_opt eq "none")
{
    $ignore_flag = "neither";
}
elsif($ignore_opt eq "both")
{
    $ignore_flag = "both";
}
else
{
    print ("INVALID -I option\n");
    print $usage;
    print $options_usage;
    exit();
}

open(ALL, ">" . $out_root . ".all");
open(FAM, ">" . $out_root . ".fam");
open(ROC, ">" . $out_root . ".roc");

print FAM "RMARK benchmark (processed with rmark_process_glbf.pl)\n";
print FAM "    module     = $rmm\n";
print FAM "    configfile = $rmk\n"; 
print FAM "    index      = $idx\n";
print FAM "    genome     = $genome_file\n";
print FAM "    glbf       = $glbf_file\n";
print FAM "    mode       = $res_opt\n\n";

print ALL "RMARK benchmark (processed with rmark_process_glbf.pl)\n";
print ALL "    module     = $rmm\n";
print ALL "    configfile = $rmk\n"; 
print ALL "    index      = $idx\n";
print ALL "    genome     = $genome_file\n";
print ALL "    glbf       = $glbf_file\n";
print ALL "    mode       = $res_opt\n\n";

#BLAST scores are E values where lower is better, this messes up
#our calculation of $highnoise, so we have a flag
$lower_better = 0;
if($rmm =~ m/blast/) { $lower_better = 1; }

if($lower_better) { $worst_score = 10; }
else {$worst_score = -1;}

@fam_roc_pos_arr = ();
@fam_roc_sc_arr = ();
@fam_roc_x_arr = ();
@fam_roc_y_arr = ();
@roc_pos_arr = ();
@roc_sc_arr = ();
@roc_x_arr = ();
@roc_y_arr = ();

# We're going to keep track of the summed MER stats across families
# these are the best possible MER scores; i.e. what we'd get if
# we had perfect E-values.
$summed_across_fam_mer_score = 0;  
$summed_across_fam_mer_fp    = 0;
$summed_across_fam_mer_fn    = 0;

# Collect information from the genome and embed files that
# tell us about where the positives are.
#
# If we are in fnt mode:
# Each nt not spanned by a positive is potentially a negative.
# Each nt that's part of a true test sequence is considered a positive
#
# If we are in nnt mode:
# Each nt not spanned by a positive is potentially a negative.
#
# %nt_HAH:      a hash of arrays of hashes will keep track of the negatives
#               and positives across all families
# %pos_nt_HAH:  a hash of arrays of hashes that will keep track of which
#               nucleotides are positives
# %fam_nt_HHAH: a hash of hashes of arrays of hashes will keep track of the
#               positives and negatives within a family.
#
# %nt_HAH:
# key 1D: name of each chromosome.
# size of array 2D: exactly 2 elements
# array 2D values : element 0 for forward direction element 1 for reverse direction.
# key 3D: nt position for key 1D orientation array 2D value
# value 3D : 'P' is part of a positive hit
#            or a number indicating the best score of any hit to that position.
#            if a given position (key) doesn't exist, this means no hit has
#            been reported to it, and it is an inferred true negative.
#
# %fam_nt_HAH is the same as nt_hah but the first dimension is keyed by the family name
#
# %pos_nt_HAH:
# key 1D: name of each chromosome.
# size of array 2D: exactly 2 elements
# array 2D values : element 0 for forward direction element 1 for reverse direction.
# key 3D: nt position for key 1D orientation array 2D value
# value 3D : family name X of positive indicating this nt position in this orientation and
#            chromosome is a positive for family X.
#

# Read the genome, just so we know num chromosomes and length of chromosomes.
%seq_hash = ();
$genome_file;
read_fasta($genome_file, \%seq_hash);
foreach $key (keys(%seq_hash))
{
    $genome_length_FR += 2*length($seq_hash{$key});
}

# Read the family list
%fam_hash = ();
open (INDEX,$idx) || die;
while (<INDEX>) 
{
    if (/^(\S+)/) {
	$curr_fam = $1;
	$fam_hash{$curr_fam} = 1;
	if(!(exists($fam_nt_HHAH{$curr_fam})))
	{
	    %{$fam_nt_HHAH{$curr_fam}} = ();
	    #print("adding for $curr_fam\n");
	    foreach $chrom_key (keys(%seq_hash))
	    {
		@{$fam_nt_HHAH{$curr_fam}{$chrom_key}} = ();
		%{$fam_nt_HHAH{$curr_fam}{$chrom_key}[0]} = ();
		%{$fam_nt_HHAH{$curr_fam}{$chrom_key}[1]} = ();
		#print("added fam_nt_HHAH{$curr_fam}{$chrom_key}\n");
	    }
	}
    }
}
close(INDEX);

#print("total non pos chars = $all_non_pos_chars\n");

# Read the embed file and determine where all the embedded sequences 
# (true positives) are in the pseudo-genome as well as which family 
# each sequence belongs to.
# We take advantage of fact that true positives cannot overlap
# based on how we constructed the pseudo-genome.
open(EBDLIST,"$embed_file") || die;
while (<EBDLIST>) {
    if (/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) { 
	$curr_fam = $1;
	$curr_tp = $2;
	$curr_chrom = $3;
	$curr_begin = $4;
	$curr_end = $5;
	$curr_orient = $6;
	#we only care about families from the fam.list file
	#print("curr fam : $curr_fam\n");
	if(exists($fam_hash{$curr_fam}))
	{
	    #print("adding positives\n");
	    $pos_chrom_begin_hash{$curr_fam}{$curr_orient}{$curr_chrom}{$curr_tp} = $curr_begin;
	    $pos_chrom_end_hash{$curr_fam}{$curr_orient}{$curr_chrom}{$curr_tp} = $curr_end;
	    $haspos{$curr_chrom} = 1;
	    $ntrue_perfam{$curr_fam}++;
	    $ntrue_perchrom{$curr_chrom}++;
	    $ntrue_all++;
	    
	    # Add the information on this positive to %nt_HHAH (if we're in NT mode)
	    for($i = $curr_begin; $i <= $curr_end; $i++)
	    {
		$pos_nt_HAH{$curr_chrom}[$curr_orient]{$i} = $curr_fam;
		$pos_chars_per_fam{$curr_fam}++;
		$all_pos_chars++;
	    }
	}
    }
}
close EBDLIST;


# Read the *.glbf file, and parse it appropriately. 
# We have a separate function for each mode.

#In these functions we fill data structures
#that its easy to sort by score with. 
#The old structure
#%fam_nt_hhah was used just to ease implementation.
#
#fam_sc_HH is a 2D hash:
# 1D key = family name
# 2D key = hit name ($family . "." . $chrom . "." . $orientation . "." position)
#    value = best score to this nt
#all_sc_H is a 1D hash:
#    key = hit name ($family . "." . $chrom . "." . $orientation . "." position)
#    value = best score to this nt
#all_is_pos_H is a 1D hash:
#    key = hit name (same as in all_sc_H and fam_sc_HH

if($fnt)
{
    parse_glbf_fnt($glbf_file, \%fam_nt_HHAH, \%pos_nt_HAH, \%fam_sc_HH, \%all_sc_H, 
		   \%all_is_pos_H, $lower_better, $ignore_flag);
}
elsif($nnt)
{
    parse_glbf_nnt($glbf_file, \%fam_nt_HHAH, \%pos_nt_HAH, \%pos_chrom_begin_hash, 
		   \%pos_chrom_end_hash, \%fam_sc_HH, \%all_sc_H, \%all_is_pos_H, 
		   $lower_better, $ignore_flag);
}
elsif($hit)
{
    parse_glbf_hit($glbf_file, \%pos_chrom_begin_hash, 
		   \%pos_chrom_end_hash, \%fam_sc_HH, \%all_sc_H, \%all_is_pos_H, 
		   $lower_better, $ignore_flag);
}

# For each family, write out the results, do this according to the order in the index file:

open (INDEX,$idx) || die;
while (<INDEX>) 
{
    if (/^(\S+)/) 
    {
	$fam = $1;
	if($fnt)
	{
	    printf FAM ("family: %-15s (%3d pos seqs covering %d nts)\n", $fam, $ntrue_perfam{$fam}, $pos_chars_per_fam{$fam});
	}
	else
	{
	    printf FAM ("family: %-15s (%3d positives)\n", $fam, $ntrue_perfam{$fam});
	}

        #sort all the family hits by score, and print them out in sorted order
	$sawnoise = 0;
	$num_above_noise = 0;
	$fam_nneg = 0;
	$fam_npos = 0;
	#keep track of negative characters for NT mode
	$all_non_pos_chars += $genome_length_FR;
	if($ignore_flag eq "neither")
	{
	    $all_non_pos_chars -= $pos_chars_per_fam{$fam};
	    $fam_non_pos_chars = $genome_length_FR - $pos_chars_per_fam{$fam};
	}
	elsif($ignore_flag eq "opposite")
	{
	    #if we're ignoring cross hits on the other strand, 
	    #then we subtract ALL the positive characters
	    $all_non_pos_chars -= $all_pos_chars;
	    $fam_non_pos_chars = $genome_length_FR - $pos_chars_per_fam{$fam};
	}
	elsif($ignore_flag eq "both")
	{
	    #if we're ignoring cross hits on either strand, 
	    #then we subtract 2 * ALL the positive characters
	    $all_non_pos_chars -= 2 * $all_pos_chars;
	}
	if($lower_better)
	{
	    @sorted_hits = sort {$fam_sc_HH{$fam}{$a} <=> $fam_sc_HH{$fam}{$b}} keys(%{$fam_sc_HH{$fam}});
	}
	else
	{
	    @sorted_hits = sort {$fam_sc_HH{$fam}{$b} <=> $fam_sc_HH{$fam}{$a}} keys(%{$fam_sc_HH{$fam}});
	}
	@fam_roc_pos_arr = ();
	@fam_roc_sc_arr = ();
	@fam_roc_x_arr = ();
	@fam_roc_y_arr = ();
	
	foreach $hit (@sorted_hits)
	{
	    if(!(exists($all_is_pos_H{$hit}))) { 
		$plus_or_neg = "-"; 
		$fam_nneg++;
		if(!($sawnoise)) {
		    $sawnoise = 1;
		    $num_above_noise = $fam_npos;
		}
	    }
	    else {$plus_or_neg = "+"; $fam_npos++; }
	    printf FAM ("%s %-40s %g\n", $plus_or_neg, $hit, $fam_sc_HH{$fam}{$hit});
	    #printf("%s %-40s %g\n", $plus_or_neg, $hit, $fam_sc_HH{$fam}{$hit});
	    push(@fam_roc_pos_arr, (exists($all_is_pos_H{$hit})));
	    push(@fam_roc_sc_arr, $fam_sc_HH{$fam}{$hit});
	}
	if(!($sawnoise)) {
	    $num_above_noise = $fam_npos;
	}
	if($fnt)
	{
	    mer_and_roc_points(\@fam_roc_sc_arr, \@fam_roc_pos_arr, \@fam_roc_x_arr, \@fam_roc_y_arr, $pos_chars_per_fam{$fam}, $fam_non_pos_chars, \$mer_thresh, \$mer_score, \$mer_fp, \$mer_fn);
	}
	elsif($nnt)
	{
	    mer_and_roc_points(\@fam_roc_sc_arr, \@fam_roc_pos_arr, \@fam_roc_x_arr, \@fam_roc_y_arr, $ntrue_perfam{$fam}, $fam_non_pos_chars, \$mer_thresh, \$mer_score, \$mer_fp, \$mer_fn);
	}
	elsif($hit)
	{
	    mer_and_roc_points(\@fam_roc_sc_arr, \@fam_roc_pos_arr, \@fam_roc_x_arr, \@fam_roc_y_arr, $ntrue_perfam{$fam}, $fam_nneg, \$mer_thresh, \$mer_score, \$mer_fp, \$mer_fn);
	}
	printf FAM ("positives:  %d | above noise: %d\n", $fam_npos, $num_above_noise);
	printf FAM ("negatives:  %d\n", $fam_nneg);
	printf FAM ("MER:        %d\n", $mer_score);
	printf FAM ("MER_fp:     %d\n", $mer_fp);
	printf FAM ("MER_fn:     %d\n", $mer_fn);
	printf FAM ("MER_thresh: %g\n\n", $mer_thresh);
	$summed_across_fam_mer_score += $mer_score;
	$summed_across_fam_mer_fp    += $mer_fp;
	$summed_across_fam_mer_fn    += $mer_fn;
    }
}
# Next sort all the scores, and print the results to the *.all file.
$sawnoise = 0;
printf FAM ("\n\nTotal: (%3d pos seqs covering %d nts)\n", $ntrue_all, $all_pos_chars);
printf ALL ("\n\nTotal: (%3d pos seqs covering %d nts)\n", $ntrue_all, $all_pos_chars);
if($lower_better)
{
    @sorted_hits = sort {$all_sc_H{$a} <=> $all_sc_H{$b}} keys(%all_sc_H);
}
else
{
    @sorted_hits = sort {$all_sc_H{$b} <=> $all_sc_H{$a}} keys(%all_sc_H);
}
@roc_pos_arr = ();
@roc_sc_arr = ();
foreach $hit (@sorted_hits)
{
    if(!(exists($all_is_pos_H{$hit}))) 
    { 
	$plus_or_neg = "-"; 
	$nneg++;
	if(!($sawnoise)) {
	    $sawnoise = 1;
	    $num_above_noise = $npos;
	}
    }
    else {$plus_or_neg = "+"; $npos++;}
    printf ALL ("%s %-40s %g\n", $plus_or_neg, $hit, $all_sc_H{$hit});
    push(@roc_pos_arr, (exists($all_is_pos_H{$hit})));
    push(@roc_sc_arr, $all_sc_H{$hit});
    #print("pushed $hit sc : $all_sc_H{$hit} to roc_sc_arr\n");
}
#now we add the final values to the roc_pos_arr and roc_sc_arr which
#denote the cutpoint if we set the threshold to the lowest possible 
#score, such that EVERY nt was considered a positive.

# now determine the ROC curve points based on all families combined (full sorted list)
@roc_x_arr = ();
@roc_y_arr = ();
if($fnt)
{
    mer_and_roc_points(\@roc_sc_arr, \@roc_pos_arr, \@roc_x_arr, \@roc_y_arr, $all_pos_chars, $all_non_pos_chars, \$mer_thresh, \$mer_score, \$mer_fp, \$mer_fn);
}
elsif($nnt)
{
    mer_and_roc_points(\@roc_sc_arr, \@roc_pos_arr, \@roc_x_arr, \@roc_y_arr, $ntrue_all, $all_non_pos_chars, \$mer_thresh, \$mer_score, \$mer_fp, \$mer_fn);
}
elsif($hit)
{
    mer_and_roc_points(\@roc_sc_arr, \@roc_pos_arr, \@roc_x_arr, \@roc_y_arr, $ntrue_all, $nneg, \$mer_thresh, \$mer_score, \$mer_fp, \$mer_fn);
}
printf ALL ("positives: %d | above noise: %d\n", $npos, $num_above_noise );
printf ALL ("negatives: %d\n\n", $nneg);
printf FAM ("positives: %d | above noise: %d\n", $npos, $num_above_noise );
printf FAM ("negatives: %d\n\n", $nneg);

printf ALL ("MER:        %d\n", $mer_score);
printf ALL ("MER_fp:     %d\n", $mer_fp);
printf ALL ("MER_fn:     %d\n", $mer_fn);
printf ALL ("MER_thresh: %g\n\n", $mer_thresh);

printf FAM ("MER:        %d\n", $mer_score);
printf FAM ("MER_fp:     %d\n", $mer_fp);
printf FAM ("MER_fn:     %d\n", $mer_fn);
printf FAM ("MER_thresh: %g\n\n", $mer_thresh);

printf FAM ("MER statistics summed across all " . scalar(keys(%fam_hash)) . " families:\n");
printf FAM ("MER    (fam sum):     %d\n", $summed_across_fam_mer_score);
printf FAM ("MER_fp (fam sum):     %d\n", $summed_across_fam_mer_fp);
printf FAM ("MER_fn (fam sum):     %d\n", $summed_across_fam_mer_fn);

for($i = 0; $i < scalar(@roc_x_arr); $i++)
{
    printf ROC ("(" . $roc_x_arr[$i] . ", " . $roc_y_arr[$i] . ")\n");
}
close(FAM);
close(ALL);
close(ROC);

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
	return 1;
    }
    return 0;
}
#################################################################

#################################################################
# subroutine : mer_and_roc_points
#
# EPN 09.14.05
#
# purpose : Generate (x,y) data points for an ROC curve x=FP fraction
#           y=TP fraction.  Given a list of hits sorted by score and 
#           information specifying each hit as positive or negative, 
#           as well number of positives.  Go through list and generate
#           a new data point for member of the list.
#          
#           For each (x,y) point 
#            o x = 1-specificity = 1-(TN/(TN+FP)) 
#            o y = sensitivity   =    TP/(TP+FN)
#            for some given cut-point
#            see: http://gim.unmc.edu/dxtests/roc1.htm
#                 for a good reference (intro to ROC curves)
#
#           (09.25.05) Also determine the MER threshold and score 
#           (minimum error rate).
#           MER thresholed : (a cutpoint in the ranked list) at which 
#           the sum of FP and FN is minimized.  The MER score is (FP+FN)
#           at the MER threshold.
#
# args : (1) $sc_arr_ref
#            reference to array with scores for each hit
#        (2) $pos_arr_ref
#            reference to array where arr[$i]=1 indicates hit $i
#            is a positive, and arr[$i]=0 indicates hit $i is
#            a negative.
#        (3) $x_arr_ref
#            array to fill with x coordinates of roc points
#        (4) $y_arr_ref
#            array to fill with y coordinates of roc points
#        (5) $npos
#            total number of positives (TP+FN)
#        (6) $nneg
#            total number of negatives (TN+FP)
#        (7) $mer_thresh_ref
#            ref to scalar that will hold the MER threshold
#        (8) $mer_score_ref
#            ref to scalar that will hold the MER score (FP+FN at
#            $mer_thresh
#        (9) $mer_fp_ref
#            ref to scalar that will be the number of fp at MER threshold
#       (10) $mer_fn_ref
#            ref to scalar that will be the number of fn at MER threshold
################################################################# 
sub mer_and_roc_points
{
    if(scalar(@_) != 10)
    {
	die "ERROR in mer_and_roc_points, exactly 10 arguments expected.\n";
    }

    my($sc_arr_ref, $pos_arr_ref, $x_arr_ref, $y_arr_ref, $npos, $nneg,
       $mer_thresh_ref, $mer_score_ref, $mer_fp_ref, $mer_fn_ref) = @_;

    $nhit = scalar(@{$sc_arr_ref});
    #print("nhit : $nhit\n");
    #set up our false negative, true negative, false positive, true positive counts
    $fn = $npos;
    $tn = $nneg;
    #printf("init fn : $npos | tn : $nneg\n");
    $fp = 0;
    $tp = 0;
    $prev_was_pos = 1;
    $prev_was_neg = 1;
    $mer_thresh = $sc_arr_ref->[0];
    $mer_score = $fn + $fp;
    $mer_fp = $fp;
    $mer_fn = $fn;
    for($i = 0; $i < ($nhit-1); $i++)
    {
	$next_is_pos = $pos_arr_ref->[($i+1)];
	$next_is_neg = (!($pos_arr_ref->[($i+1)]));
    	
	if($pos_arr_ref->[$i])
	{
	    #positive
	    $tp++;
	    $fn--;
	    if($prev_was_neg || $next_is_neg)
	    {
		$add_point = 1;
	    }
	    else
	    {
		$add_point = 0;
	    }
	    $prev_was_pos = 1;
	    $prev_was_neg = 0;
	}
	elsif(!($pos_arr_ref->[$i]))
	{
	    #negative
	    $fp++;
	    $tn--;
	    if($prev_was_pos || $next_is_pos)
	    {
		$add_point = 1;
	    }
	    else
	    {
		$add_point = 0;
	    }
	    $prev_was_pos = 0;
	    $prev_was_neg = 1;
	}
	if($add_point)
	{
	    if(($tn+$fp) > 0)
	    {
		push(@{$x_arr_ref}, (1-($tn/($tn+$fp))));
	    }
	    else
	    {
		push(@{$x_arr_ref}, 1);
	    }
	    if(($tp+$fn) > 0)
	    {
		push(@{$y_arr_ref}, ($tp/($tp+$fn)));
	    }
	    else
	    {
		push(@{$y_arr_ref}, 0);
	    }
	}
	#check if we have a new MER
	if(($mer_score) > ($fp+$fn))
	{
	    $mer_score = $fp+$fn;
	    $mer_thresh = $sc_arr_ref->[($i+1)];
	    $mer_fp = $fp;
	    $mer_fn = $fn;
	}
    }
    #now add the final point
    if($pos_arr_ref->[($nhit-1)])
    {
	#positive
	$tp++;
	$fn--;
    }
    elsif(!($pos_arr_ref->[($nhit-1)]))
    {
	#negative
	$fp++;
	$tn--;
    }
    if(($tn+$fp) > 0)
    {
	push(@{$x_arr_ref}, (1-($tn/($tn+$fp))));
    }
    else
    {
	push(@{$x_arr_ref}, 1);
    }
    if(($tp+$fn) > 0)
    {
	push(@{$y_arr_ref}, ($tp/($tp+$fn)));
    }
    else
    {
	push(@{$y_arr_ref}, 0);
    }
    #we could have an MER_threshold that is 'below' the final
    #score on the ranked list, we'll handle this specially
    if(($mer_score) > ($fp+$fn))
    {
	$mer_score = $fp+$fn;
	#assumes HIGHER score is better
	$mer_thresh = "<" . $sc_arr_ref->[($i)];
	$mer_fp = $fp;
	$mer_fn = $fn;
    }
    $$mer_score_ref = $mer_score;
    $$mer_thresh_ref = $mer_thresh;
    $$mer_fp_ref = $mer_fp;
    $$mer_fn_ref = $mer_fn;
}
################################################################# 

#################################################################
# from       : M_seq.pm
# subroutine : read_fasta
# sub class  : crw and sequence
# 
# EPN 03.08.05
#
# purpose : Open, read, and store the information in a given
#           .fa (fasta format) file.
#
# args : (1) $in_file
#            name of .fa file in current directory
#        (2) $seq_hash_ref
#            reference to the hash that will contain the sequence
#            information.  Fasta description line used as key for
#            each sequence, sequence is value.
################################################################# 
sub read_fasta
{
    my($in_file, $seq_hash_ref) = @_;
    open(IN, $in_file) || die "ERROR in main::read_fasta() could not open $in_file";
    
    #chomp up beginning blank lines
    $line = <IN>;
    while(!($line =~ m/^>/))
    {
	$line = <IN>;
    }

    chomp $line;
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
	$seq_name = $line;
	$seq_name =~ s/^>//;
    }
    #trim_keys_in_hash($seq_hash_ref, DEFAULT_MAX_SEQ_HEADER_LENGTH);
}

#################################################################

#################################################################
# 
# subroutine : parse_glbf_fnt
# sub class  : sequence
# 
# EPN 10.11.05
#
# purpose : Open, read, and store the information in a given
#           .glbf file using the fnt (FULL NT) resolution mode.
#
# args : (1) $in_file
#            name of .glbf file in current directory
#        (2) $fam_nt_HHAHR
#            reference to a hash of hashes of hashes of arrays of hashes
#            key 1D: family name
#            key 2D: name of each chromosome.
#            size of array 3D: exactly 2 elements
#            array 3D values : element 0 for forward direction element 1 for reverse direction.
#            key 4D: nt position for key 1D orientation array 2D value
#            value 4D :  the best score of any hit to that position.
#                        if a given position (key) doesn't exist, this means no hit has
#                        been reported to it, and it is an inferred true negative.
#        (3) $pos_nt_HAHR
#            reference to a HAH
#            key 1D: name of each chromosome.
#            size of array 2D: exactly 2 elements
#            array 2D values : element 0 for forward direction element 1 for reverse direction.
#            key 3D: nt position for key 1D orientation array 2D value
#            value 3D : family name X of positive indicating this nt position in this orientation and
#            chromosome is a positive for family X.
#        (4) $fam_sc_HHR
#            refererence to a HH (2D hash) TO BE FILLED IN THIS SUB
#            key 1D: name of family
#            key 2D: name of hit ($fam "." $chrom "." $orient "." $position) 
#            value : score of hit
#        (5) $all_sc_HR
#            reference to a hash TO BE FILLED IN THIS SUB
#            key   : name of hit ($fam "." $chrom "." $orient "." $position) 
#            value : score of hit
#        (6) $all_is_pos_HR
#            reference to a hash TO BE FILLED IN THIS SUB
#            key   : name of A POSITIVE hit
#            value : 1 (always)
#        (7) $lower_better
#            1 if a lower score is better, 0 if higher is better.
#        (8) $ignore_flag
#            "opposite" to ignore any 'negative' hit from family X to 
#            a true test sequence in family Y, where X != Y on opposite strand only.
#            "both" to ignore any 'negative' hit from family X to 
#            a true test sequence in family Y, where X != Y on both strands.
#            "neither" to not ignore any cross hits.
#            
################################################################# 
sub parse_glbf_fnt
{
    if(scalar(@_) != 8)
    {
	die "ERROR in main::parse_glbf_fnt(), expecting 8 arguments\n";
    }
    my($in_file, $fam_nt_HHAHR, $pos_nt_HAHR, $fam_sc_HHR, $all_sc_HR, $all_is_pos_HR,
       $lower_better, $ignore_flag) = @_;
    open(IN, $in_file) || die "ERROR in main::parse_glbf_fnt() could not open $in_file";

    while($line = <IN>)
    {
	if ($line =~ /\>(\S+).*/)
	{
	    $fam = $1;
	}
	elsif ($line =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$/)
	{
	    $chrom = $1;
	    $score = $2;
	    $begin = $3;
	    $end = $4;
	    $orient = $5;
	    
	    #its possible that infernal reported the hit as a reverse so
	    #$begin > $end, switch these two in this case.
	    if($begin > $end) { $temp = $end; $end = $begin; $begin = $temp; }
	    
	    # Fill in our %nt_HAH and %fam_nt_HHAH with the information
	    # for each nucleotide in the current hit.
	    #
	    # pseudo code:
	    # for each nt in the hit
	    #   check to see if a previous hit spanned it.  
		#   If yes, we check to see if the old hit's score was better.
	    #      If yes, we do nothing.
	    #      If no, replace the score with our new score.
	    #   If no, this is the first hit to span that nt.
	    #      Set score with our new score.
	    #
	    # we'll deal with determining if the hits are positives
	    # or negatives later.
			
	    for ($i = $begin; $i <= $end; $i++)
	    {
		#print("\ti : $i | chrom $chrom\n");
		
		if(!(exists($fam_nt_HHAHR->{$fam}{$chrom})))
		{
		    die "ERROR in rm_process_glbf.pl::main we've got a hit to a chromosome that we didn't even know existed!\n";
		}
		if(exists($fam_nt_HHAHR->{$fam}{$chrom}[$orient]{$i}))
		{
		    #we've already seen a hit to this nt
		    $old_sc = $fam_nt_HHAHR->{$fam}{$chrom}[$orient]{$i};
		    if((($lower_better) && 
			($old_sc < $score))
		       || 
		       ((!($lower_better)) && 
			($old_sc > $score)))
		    {
			#old hit's score is better - do nothing 
		    }
		    else
		    {
			#new score is better, replace old score
			$fam_nt_HHAHR->{$fam}{$chrom}[$orient]{$i} = $score;
		    }
		}
		else
		{
		    #first hit to this nt - store score
		    $fam_nt_HHAHR->{$fam}{$chrom}[$orient]{$i} = $score;
		}
	    }
	}
    }
    #now name all the hits and add them all by name to 
    #the fam_pos_neg_sc_HH and all_pos_neg_sc_H
    #, this will be a data structure
    #that its easy to sort by score with. The old structure
    #%fam_nt_hhah was used just to ease implementation.
    #
    #fam_sc_HH is a 2D hash:
    # 1D key = family name
    # 2D key = hit name ($family . "." . $chrom . "." . $orientation . "." position)
    #    value = best score to this nt
    #all_sc_H is a 1D hash:
    #    key = hit name ($family . "." . $chrom . "." . $orientation . "." position)
    #    value = best score to this nt
    #all_is_pos_H is a 1D hash:
    #    key = hit name (same as in all_sc_H and fam_sc_HH
    foreach $fam (keys(%{$fam_nt_HHAHR}))
    {
	foreach $chrom (keys(%{$fam_nt_HHAHR->{$fam}}))
	{
	    for($o = 0; $o < scalar(@{$fam_nt_HHAHR->{$fam}{$chrom}}); $o++)
	    {
		#$o = 0 or 1 (orientation)
		foreach $posn (keys(%{$fam_nt_HHAHR->{$fam}{$chrom}[$o]}))
		{
		    #print("posn : $posn\n");
		    $score = $fam_nt_HHAHR->{$fam}{$chrom}[$o]{$posn};
		    $new_key = $fam . "." . $chrom . "." . $o . "." . "$posn";
		    #Add to family specific hash
		    $fam_sc_HHR->{$fam}{$new_key} = $score;
		    #add to master hash
		    $all_sc_HR->{$new_key} = $score;
		    #printf("added hit to fam $fam : $new_key\n");
		    #now if its a positive, add it to all_is_pos_H
		    if(exists($pos_nt_HAHR->{$chrom}[$o]{$posn}) &&
		       $pos_nt_HAHR->{$chrom}[$o]{$posn} eq $fam)
		    {
			$all_is_pos_HR->{$new_key} = 1;
		    }
		    elsif(exists($pos_nt_HAHR->{$chrom}[$o]{$posn}) &&
		       $pos_nt_HAHR->{$chrom}[$o]{$posn} ne $fam)
		    {
			#new added 10.17.05
			#IF WE'RE IGNORING CROSS-FAMILY HITS (hits from family X to family Y)
			#we want to ignore this hit!
			if(($ignore_flag eq "opposite" || $ignore_flag eq "both"))
			{
			    delete($fam_sc_HHR->{$fam}{$new_key});
			    delete($all_sc_HR->{$new_key});
			}
		    }	
		    #new added 10.28.05
		    #If we're ignoring cross-family hits on either strand, 
		    #and this hit overlaps with any true test sequence
		    #on the other strand, we ignore this hit.
		    if($ignore_flag eq "both")
		    {
			if($o == 0) { $other_o = 1;}
			else { $other_o = 0;}
			if(exists($pos_nt_HAHR->{$chrom}[$other_o]{$posn}))
			{
			    delete($fam_sc_HHR->{$fam}{$new_key});
			    delete($all_sc_HR->{$new_key});
			}
		    }
		}
	    }
	}
    }
}
#################################################################

#################################################################
# 
# subroutine : parse_glbf_nnt
# sub class  : sequence
# 
# EPN 10.11.05
#
# purpose : Open, read, and store the information in a given
#           .glbf file using the nnt resolution mode.
#
# args : (1) $in_file
#            name of .glbf file in current directory
#        (2) $fam_nt_HHAHR
#            reference to a hash of hashes of hashes of arrays of hashes
#            key 1D: family name
#            key 2D: name of each chromosome.
#            size of array 3D: exactly 2 elements
#            array 3D values : element 0 for forward direction element 1 for reverse direction.
#            key 4D: nt position for key 1D orientation array 2D value
#            value 4D :  the best score of any hit to that position.
#                        if a given position (key) doesn't exist, this means no hit has
#                        been reported to it, and it is an inferred true negative.
#        (3) $pos_nt_HAHR
#            reference to a HAH
#            key 1D: name of each chromosome.
#            size of array 2D: exactly 2 elements
#            array 2D values : element 0 for forward direction element 1 for reverse direction.
#            key 3D: nt position for key 1D orientation array 2D value
#            value 3D : family name X of positive indicating this nt position in this orientation and
#            chromosome is a positive for family X.
#        (4) $pos_chrom_begin_HHHHR
#            reference to a 4D hash
#            key 1D: name of a family
#            key 2D: name of orientation
#            key 3D: name of chromosome
#            key 4D: name of true test sequence (name from Rfam alignment)
#            value : beginning position of true test sequence
#        (5) $pos_chrom_end_HHHHR
#            reference to a 4D hash
#            key 1D: name of a family
#            key 2D: name of orientation
#            key 3D: name of chromosome
#            key 4D: name of true test sequence (name from Rfam alignment)
#            value : end position of true test sequence
#        (6) $fam_sc_HHR
#            refererence to a HH (2D hash) TO BE FILLED IN THIS SUB
#            key 1D: name of family
#            key 2D: name of hit ($fam "." $chrom "." $orient "." $position) 
#            value : score of hit
#        (7) $all_sc_HR
#            reference to a hash TO BE FILLED IN THIS SUB
#            key   : name of hit ($fam "." $chrom "." $orient "." $position) 
#            value : score of hit
#        (8) $all_is_pos_HR+
#            reference to a hash TO BE FILLED IN THIS SUB
#            key   : name of A POSITIVE hit
#            value : 1 (always)
#        (9) $lower_better
#            1 if a lower score is better, 0 if higher is better.
#       (10) $ignore_flag
#            "opposite" to ignore any 'negative' hit from family X to 
#            a true test sequence in family Y, where X != Y on correct strand only.
#            "both" to ignore any 'negative' hit from family X to 
#            a true test sequence in family Y, where X != Y on both strands.
#            "neither" to not ignore any cross hits.
#
################################################################# 
sub parse_glbf_nnt
{
    if(scalar(@_) != 10)
    {
	die "ERROR in main::parse_glbf_nnt(), expecting 10 arguments\n";
    }
    my($in_file, $fam_nt_HHAHR, $pos_nt_HAHR, $pos_chrom_begin_HHHHR, $pos_chrom_end_HHHHR, $fam_sc_HHR, $all_sc_HR, $all_is_pos_HR, $lower_better, $ignore_flag) = @_;
    open(IN, $in_file) || die "ERROR in main::parse_glbf_nnt() could not open $in_file";

    $min_overlap_fract = 0.5;
    while($line = <IN>)
    {
	if ($line =~ /\>(\S+).*/)
	{
	    $fam = $1;
	}
	elsif ($line =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$/)
	{
	    $chrom = $1;
	    $score = $2;
	    $begin = $3;
	    $end = $4;
	    $orient = $5;
	    #printf("\n\nexamining hit : chrom : $chrom | score : $score | begin : $begin | end : $end\n");
	    #its possible that infernal reported the hit as a reverse so
	    #$begin > $end, switch these two in this case.
	    if($begin > $end) { $temp = $end; $end = $begin; $begin = $temp; }
	    #is this a positive hit?
	    $found_pos = 0;
	    foreach $pos_key (keys(%{$pos_chrom_begin_HHHHR->{$fam}{$orient}{$chrom}}))
	    {
		#print("pos : $pos_key | begin : $pos_chrom_begin_hash{$fam}{$chrom}{$pos_key} | end : $pos_chrom_end_hash{$fam}{$chrom}{$pos_key}\n");
		#for the overlap call its important $begin and $end of the
		#hit go first, requiring only that $min_overlap_fract
		#of the hit must be inside the positive region.
		if(overlap($begin,
			   $end,
			   $pos_chrom_begin_HHHHR->{$fam}{$orient}{$chrom}{$pos_key},
			   $pos_chrom_end_HHHHR->{$fam}{$orient}{$chrom}{$pos_key},
			   $min_overlap_fract))
		{
		    if($found_pos) #this means our hit overlaps two positives
			#not sure what to do here - abort!
		    {
			die "ERROR in parse_glbf_nnt hit to $chrom : $begin-$end overlaps > 1 positive.\n";
		    }
		    $found_pos = 1;
		    $pos_id = $pos_key;
		}
	    }
	    if($found_pos) 
	    {
		#printf("found a positive\n");
		#the current hit overlaps with a positive
		if(exists($found_pos_hash{$pos_id})) 
		{
		    #a different hit overlaps with the same positive
		    if((($lower_better) && 
			($fam_sc_HHR->{$fam}{$pos_id} < $score))
		       || 
		       ((!($lower_better)) && 
			($fam_sc_HHR->{$fam}{$pos_id} > $score)))
		    {
			#the old hit that overlaps with this positive is better!
			#do nothing
			#printf("old hit : $pos_id sc : $fam_sc_HHR->{$fam}{$pos_id} is better\n");
		    }
		    else
		    {
			#new hit has better score, modify score of best hit.
			#printf("new hit sc $score better than old hit $pos_id sc : $fam_sc_HHR->{$fam}{$pos_id}\n");
			$fam_sc_HHR->{$fam}{$pos_id} = $score;
			$all_sc_HR->{$pos_id} = $score;
			$all_is_pos_HR->{$pos_id} = $score;
		    }
		}
		else
		{
		    #there currently is no hit to this positive hit
		    #add the new one
		    $fam_sc_HHR->{$fam}{$pos_id} = $score;
		    $all_sc_HR->{$pos_id} = $score;
		    $all_is_pos_HR->{$pos_id} = $score;
		    $found_pos_hash{$pos_id} = 1;
		    #printf("no current hit to $pos_id\n");
		}
	    }
	    else
	    {
		#this hit does not overlap with a positive, its a negative
		#new added 10.17.05
		#IF WE'RE IGNORING CROSS-FAMILY HITS (hits from family X to family Y)
		#we IGNORE THIS HIT!
		$ignore_this_seq = 0;
		if($ignore_flag eq "opposite")
		{
		    foreach $other_fam(keys(%{$pos_chrom_begin_HHHHR}))
		    {
			if($other_fam ne $fam)
			{
			    foreach $pos_key (keys(%{$pos_chrom_begin_HHHHR->{$other_fam}{$orient}{$chrom}}))
			    {
				if(overlap($begin,
					   $end,
					   $pos_chrom_begin_HHHHR->{$other_fam}{$orient}{$chrom}{$pos_key},
					   $pos_chrom_end_HHHHR->{$other_fam}{$orient}{$chrom}{$pos_key},
					   $min_overlap_fract))
				{
				    $ignore_this_seq = 1;
				}
			    }
			}
		    }	
		}

		#this hit does not overlap with a positive, its a negative
		#new added 10.28.05
		#IF WE'RE IGNORING CROSS-FAMILY HITS (hits from family X to family Y)
		#ON EITHER STRAND we IGNORE THIS HIT!
		elsif($ignore_flag eq "both")
		{
		    foreach $other_fam(keys(%{$pos_chrom_begin_HHHHR}))
		    {
			foreach $orient_key (0,1)
			{
			    foreach $pos_key (keys(%{$pos_chrom_begin_HHHHR->{$other_fam}{$orient_key}{$chrom}}))
			    {
				if(overlap($begin,
					   $end,
					   $pos_chrom_begin_HHHHR->{$other_fam}{$orient_key}{$chrom}{$pos_key},
					   $pos_chrom_end_HHHHR->{$other_fam}{$orient_key}{$chrom}{$pos_key},
					   $min_overlap_fract))
				{
				    $ignore_this_seq = 1;
				}
			    }
			}	
		    }
		}

		#now if $ignore_this_seq is 1, we know we have a cross-hit, so we skip the rest
		if(!($ignore_this_seq))
		{	
		    #print("not a positive but a negative\n");
		    # We've got a negative.
		    # Pseudo-code of how to deal with negatives:
		    # for each nucleotide that this negative hit spans
		    #   check to see if a previous negative hit spanned it.  
		    #   If yes, we check to see if the old hit's score was better.
		    #   {
		    #      If yes, we do nothing.
		    #      If no, replace the score with our new score.
		    #   }
		    #   If no, this is the first negative to span that nt.
		    #      Set score with our new score.
		    #
		    # There is a 'special' case to handle.
		    # If the negative hit bleeds into a positive region
		    # but the hit doesn't qualify as a positive by our
		    # overlap criterion (if it did we wouldn't have entered
		    # this else{ ), then we ignore those nucleotides that bled
		    # into the positive (mainly because I'm not sure how to deal
		    # with them, they're not really negatives...)
		    #
		    
		    for ($i = $begin; $i <= $end; $i++)
		    {
			#print("\ti : $i | chrom $chrom\n");
			
			if(!(exists($fam_nt_HHAHR->{$fam}{$chrom})))
			{
			    die "ERROR in parse_glbf_nnt() we've got a hit to a chromosome that we didn't even know existed!\n";
			}
			if(exists($fam_nt_HHAHR->{$fam}{$chrom}[$orient]{$i}))
			{
			    #either we've already seen a negative hit to this nt
			    #or its spanned by a positive
			    $old_neg_sc = $fam_nt_HHAHR->{$fam}{$chrom}[$orient]{$i};
			    if($pos_nt_HAHR->{$chrom}[$orient]{$i} ne $fam)
			    {
				#its not a positive for the current family
				#it must be a negative nt we've already seen a hit to
				if((($lower_better) && 
				    ($old_neg_sc < $score))
				   || 
				   ((!($lower_better)) && 
				    ($old_neg_sc > $score)))
				{
				    #old hit's score is better
				    #do nothing
				}
				else
				{
				    #new score is better, replace old score
				    $fam_nt_HHAHR->{$fam}{$chrom}[$orient]{$i} = $score;
				}
			    }
			    else
			    {
				#curr nt is spanned by a positive, skip.
			    }
			}
			else
			{
			    #no previous negative spanned this nt and its ! a positive
			    $fam_nt_HHAHR->{$fam}{$chrom}[$orient]{$i} = $score;
			}
		    }
		}
	    }
	}
    }
    #we've already added the positive hits to %all_sc_H and %fam_sc_H
    #now name all the NEGATIVE hits and add them all by name to 
    #the fam_pos_neg_sc_HH and all_pos_neg_sc_H
    #, this will be a data structure
    #that its easy to sort by score with. The old structure
    #%fam_nt_hhah was used just to ease implementation.
    #
    #fam_sc_HH is a 2D hash:
    # 1D key = family name
    # 2D key = hit name ($family . "." . $chrom . "." . $orientation . "." position)
    #    value = best score to this nt
    #all_sc_H is a 1D hash:
    #    key = hit name ($family . "." . $chrom . "." . $orientation . "." position)
    #    value = best score to this nt
    #all_is_pos_H is a 1D hash:
    #    key = hit name NAME OF TRUE TEST SEQUENCE FROM RFAM
    foreach $fam (keys(%{$fam_nt_HHAHR}))
    {
	foreach $chrom (keys(%{$fam_nt_HHAHR->{$fam}}))
	{
	    for($o = 0; $o < scalar(@{$fam_nt_HHAHR->{$fam}{$chrom}}); $o++)
	    {
		#$o = 0 or 1 (orientation)
		foreach $posn (keys(%{$fam_nt_HHAHR->{$fam}{$chrom}[$o]}))
		{
		    #print("posn : $posn\n");
		    $score = $fam_nt_HHAHR->{$fam}{$chrom}[$o]{$posn};
		    $new_key = $fam . "." . $chrom . "." . $o . "." . "$posn";
		    #Add to family specific hash
		    $fam_sc_HHR->{$fam}{$new_key} = $score;
		    #add to master hash
		    $all_sc_HR->{$new_key} = $score;
		    #printf("added hit to fam $fam : $new_key\n");
		}
	    }
	}
    }
}
#################################################################

#################################################################
# 
# subroutine : parse_glbf_hit
# sub class  : sequence
# 
# EPN 10.11.05
#
# purpose : Open, read, and store the information in a given
#           .glbf file using the hit resolution mode.
#
# args : (1) $in_file
#            name of .glbf file in current directory
#        (2) $pos_chrom_begin_HHHHR
#            reference to a 4D hash
#            key 1D: name of a family
#            key 2D: name of orientation
#            key 3D: name of chromosome
#            key 4D: name of true test sequence (name from Rfam alignment)
#            value : beginning position of true test sequence
#        (3) $pos_chrom_end_HHHHR
#            reference to a 4D hash
#            key 1D: name of a family
#            key 2D: name of orientation
#            key 3D: name of chromosome
#            key 4D: name of true test sequence (name from Rfam alignment)
#            value : end position of true test sequence
#        (4) $fam_sc_HHR
#            refererence to a HH (2D hash) TO BE FILLED IN THIS SUB
#            key 1D: name of family
#            key 2D: name of hit ($fam "." $chrom "." $orient "." $position) 
#            value : score of hit
#        (5) $all_sc_HR
#            reference to a hash TO BE FILLED IN THIS SUB
#            key   : name of hit ($fam "." $chrom "." $orient "." $position) 
#            value : score of hit
#        (6) $all_is_pos_HR
#            reference to a hash TO BE FILLED IN THIS SUB
#            key   : name of A POSITIVE hit
#            value : 1 (always)
#        (7) $lower_better
#            1 if a lower score is better, 0 if higher is better.
#        (8) $ignore_flag
#            "opposite" to ignore any 'negative' hit from family X to 
#            a true test sequence in family Y, where X != Y on correct strand only.
#            "both" to ignore any 'negative' hit from family X to 
#            a true test sequence in family Y, where X != Y on both strands.
#            "neither" to not ignore any cross hits.
#
################################################################# 
sub parse_glbf_hit
{
    if(scalar(@_) != 8)
    {
	die "ERROR in main::parse_glbf_hit(), expecting 8 arguments\n";
    }
    my($in_file, $pos_chrom_begin_HHHHR, $pos_chrom_end_HHHHR, $fam_sc_HHR, $all_sc_HR, 
       $all_is_pos_HR, $lower_better, $ignore_flag) = @_;
    open(IN, $in_file) || die "ERROR in main::parse_glbf_hit() could not open $in_file";

    $min_overlap_fract = 0.5;
    while($line = <IN>)
    {
	if ($line =~ /\>(\S+).*/)
	{
	    $fam = $1;
	}
	elsif ($line =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$/)
	{
	    $chrom = $1;
	    $score = $2;
	    $begin = $3;
	    $end = $4;
	    $orient = $5;
	    #printf("\n\nexamining hit : chrom : $chrom | score : $score | begin : $begin | end : $end\n");
	    #its possible that infernal reported the hit as a reverse so
	    #$begin > $end, switch these two in this case.
	    if($begin > $end) { $temp = $end; $end = $begin; $begin = $temp; }
	    #is this a positive hit?
	    $found_pos = 0;
	    foreach $pos_key (keys(%{$pos_chrom_begin_HHHHR->{$fam}{$orient}{$chrom}}))
	    {
		#print("pos : $pos_key | begin : $pos_chrom_begin_hash{$fam}{$chrom}{$pos_key} | end : $pos_chrom_end_hash{$fam}{$chrom}{$pos_key}\n");
		#for the overlap call its important $begin and $end of the
		#hit go first, requiring only that $min_overlap_fract
		#of the hit must be inside the positive region.
		if(overlap($begin,
			   $end,
			   $pos_chrom_begin_HHHHR->{$fam}{$orient}{$chrom}{$pos_key},
			   $pos_chrom_end_HHHHR->{$fam}{$orient}{$chrom}{$pos_key},
			   $min_overlap_fract))
		{
		    if($found_pos) #this means our hit overlaps two positives
			#not sure what to do here - abort!
		    {
			die "ERROR in parse_glbf_hit hit to $chrom : $begin-$end overlaps > 1 positive.\n";
		    }
		    $found_pos = 1;
		    $pos_id = $pos_key;
		}
	    }
	    if($found_pos) 
	    {
		#printf("found a positive\n");
		#the current hit overlaps with a positive
		if(exists($found_pos_hash{$pos_id})) 
		{
		    #a different hit overlaps with the same positive
		    if((($lower_better) && 
			($fam_sc_HHR->{$fam}{$pos_id} < $score))
		       || 
		       ((!($lower_better)) && 
			($fam_sc_HHR->{$fam}{$pos_id} > $score)))
		    {
			#the old hit that overlaps with this positive is better!
			#do nothing
			#printf("old hit : $pos_id sc : $fam_sc_HHR->{$fam}{$pos_id} is better\n");
		    }
		    else
		    {
			#new hit has better score, modify score of best hit.
			#printf("new hit sc $score better than old hit $pos_id sc : $fam_sc_HHR->{$fam}{$pos_id}\n");
			$fam_sc_HHR->{$fam}{$pos_id} = $score;
			$all_sc_HR->{$pos_id} = $score;
			$all_is_pos_HR->{$pos_id} = $score;
		    }
		}
		else
		{
		    #there currently is no hit to this positive hit
		    #add the new one
		    $fam_sc_HHR->{$fam}{$pos_id} = $score;
		    $all_sc_HR->{$pos_id} = $score;
		    $all_is_pos_HR->{$pos_id} = $score;
		    $found_pos_hash{$pos_id} = 1;
		    #printf("no current hit to $pos_id\n");
		}
	    }
	    else
	    {
		#IF HIT MODE:
		@neg_id_arr = ();
		#this hit does not overlap with a positive, its a negative
		#IF WE'RE IGNORING CROSS-FAMILY HITS (hits from family X to family Y)
		#we IGNORE THIS HIT!
		$ignore_this_seq = 0;
		if($ignore_flag eq "opposite")
		{
		    foreach $other_fam(keys(%{$pos_chrom_begin_HHHHR}))
		    {
			if($other_fam ne $fam)
			{
			    foreach $pos_key (keys(%{$pos_chrom_begin_HHHHR->{$other_fam}{$orient}{$chrom}}))
			    {
				if(overlap($begin,
					   $end,
					   $pos_chrom_begin_HHHHR->{$other_fam}{$orient}{$chrom}{$pos_key},
					   $pos_chrom_end_HHHHR->{$other_fam}{$orient}{$chrom}{$pos_key},
					   $min_overlap_fract))
				{
				    $ignore_this_seq = 1;
				}
			    }
			}
		    }	
		}
		#this hit does not overlap with a positive, its a negative
		#new added 10.28.05
		#IF WE'RE IGNORING CROSS-FAMILY HITS (hits from family X to family Y)
		#ON EITHER STRAND we IGNORE THIS HIT!
		elsif($ignore_flag eq "both")
		{
		    foreach $other_fam(keys(%{$pos_chrom_begin_HHHHR}))
		    {
			foreach $orient_key (0,1)
			{
			    foreach $pos_key (keys(%{$pos_chrom_begin_HHHHR->{$other_fam}{$orient_key}{$chrom}}))
			    {
				if(overlap($begin,
					   $end,
					   $pos_chrom_begin_HHHHR->{$other_fam}{$orient_key}{$chrom}{$pos_key},
					   $pos_chrom_end_HHHHR->{$other_fam}{$orient_key}{$chrom}{$pos_key},
					   $min_overlap_fract))
				{
				    $ignore_this_seq = 1;
				}
			    }
			}	
		    }
		}


		#now if $ignore_this_seq is 1, we know we have a cross-hit, so we skip the rest
		if(!($ignore_this_seq))
		{
		    #we check to see if we already have a negative that overlaps
		    #with this hit significantly (> $min_overlap_fract).
		    $found_neg = 0;
		    foreach $neg_key (keys(%{$neg_chrom_HHHH{$fam}{$orient}{$chrom}}))
		    {
			#for the overlap call its important $begin and $end of the
			#hit go first, requiring only that $min_overlap_fract
			#of the hit must be inside the positive region.
			if(overlap($begin,
				   $end,
				   $neg_begin_hash{$fam}{$neg_key},
				   $neg_end_hash{$fam}{$neg_key},
				   $min_overlap_fract))
			{
			    #printf("found overlap with $neg_key\n");
			    #09.19.05 decided to allow multiple hits to same
			    #negative, either we delete them both or keep them both
			    $found_neg = 1;
			    push(@neg_id_arr, $neg_key);
			    #$neg_id = $neg_key;
			}
		    }
		    if($found_neg) 
		    {
			$keep_old = 0;
			foreach $neg_id (@neg_id_arr)
			{
			    #printf("curr hit overlaps with another negative we've already seen\n");
			    #the current hit is negative and overlaps with another negative
			    #we've already seen
			    if((($lower_better) && 
				($fam_sc_HHR->{$fam}{$neg_id} < $score))
			       || 
			       ((!($lower_better)) && 
				($fam_sc_HHR->{$fam}{$neg_id} > $score)))
			    {
				#printf("old hit $neg_id sc : $fam_sc_HHR->{$fam}{$neg_id} is better score than $score: do nothing\n");
				#the old hit that overlaps with this negative is better!
				$keep_old = 1;
				#do nothing
			    }
			    else
			    {
				#printf("old hit $neg_id sc : $fam_sc_HHR->{$fam}{$neg_id} is worse than sc : $score replace it\n");
			    }
			}
			if(!($keep_old))
			{
			    #new hit has better score, delete old hit(s), add new hit
			    foreach $neg_id (@neg_id_arr)
			    {
				delete($fam_sc_HHR->{$fam}{$neg_id});
				delete($all_sc_HR->{$neg_id});
				delete($neg_begin_hash{$fam}{$neg_id});
				delete($neg_end_hash{$fam}{$neg_id});
				delete($neg_chrom_HHHH{$fam}{$orient}{$chrom}{$neg_id});
			    }
			    $hit_name = $fam . "." . $chrom . "." . $orient . "." . $begin . "-" . $end;
			    $all_sc_HR->{$hit_name} = $score;
			    $fam_sc_HHR->{$fam}{$hit_name} = $score;
			    $neg_begin_hash{$fam}{$hit_name} = $begin;
			    $neg_end_hash{$fam}{$hit_name} = $end;
			    $neg_chrom_HHHH{$fam}{$orient}{$chrom}{$hit_name} = 1;
			}
		    }
		    else
		    {
			#printf("new negative from $begin to $end on $chrom\n");
			#there currently is no other negative overlapping with this negative
			#add the new one
			$hit_name = $fam . "." . $chrom . "." . $orient . "." . $begin . "-" . $end;
			$fam_sc_HHR->{$fam}{$hit_name} = $score;
			$all_sc_HR->{$hit_name} = $score;
			$neg_begin_hash{$fam}{$hit_name} = $begin;
			$neg_end_hash{$fam}{$hit_name} = $end;
			$neg_chrom_HHHH{$fam}{$orient}{$chrom}{$hit_name} = 1;
		    }
		}
	    }
	}
    }
}
