#!/usr/local/bin/perl
#
# Eric Nawrocki 10.09.05
# rmark_clusterfy.pl
#
# Prepare a RMARK benchmark for running on the cluster.
# In non-MPI mode: split a RMARK run for X families in Y 
# chromosomes into X*Y separate RMARK runs, one for each 
# family/chromosome pair.
# In MPI mode (currently only works with Infernal): 
# Split a RMARK run for X families in however many chrom-
# osomes into X separate 'mpirun' calls followed by
# a bunch of infernal2glbf calls.
#
# Usage: perl rmark_clusterfy.pl
#             <.rmm file name>
#             <.rmk file name>
#             <seq dir with *.ali *.idx *.test and *.raw files (and possibly *ali.cm files)>
#             <index file with fam names; provide path>
#             <genome root X, X.fa, X.ebd, and X.chrlist must be in seq dir>
#             <output file root>
# Options:
#        -E <x> : use E-values [default], set max E-val to keep as <x> [df: 2]
#        -B <x> : use bit scores, set min score to keep as <x>
#        -M <x> : pass -M <x> to .rmm module use, to use CM files in dir <x> instead of building them as needed
# Example:  perl rmark_clusterfy.pl infernal.rmm inf-72.rmk rmark-test/ rmark-test 
#                                   inf-72
#
# The example run above will create the following:
#     - a inf-72_rmark-test_out_dir directory with all the files needed to run
#       rmark.pl copied to it (except seq files which stay in <seq dir>)
#     - the inf-72.com file (in inf-72_rmark-test_out_dir/) a shell script 
#       which will submit the benchmark jobs to the cluster when executed
#     - the inf-72_pp.script (in inf-72_rmark-test_out_dir/) a shell script 
#       which will post-process and combine the results from all the jobs
#       to be run after all jobs finish running.
#
# General Strategy (non-MPI mode)
# (A) Create a script to execute on a cluster that will execute
#     rmark.pl many times - namely one job for each family/chromosome
#     pair in the full benchmark (SEE 'IMPORTANT 1:' below).
# (B) Copy all required files into a new directory where the jobs 
#     will be run (SEE 'IMPORTANT 2:' below).
# (C) Create a shell script to post-process the results, that will run 
#     after all the jobs have finished running.
#
# General Strategy (MPI mode)
# (A) Create a script to submit X jobs (one for each family) to the 
#     cluster to run cmsearch in MPI mode. 
# (B) Copy all required files into a new directory where the jobs 
#     will be run (SEE 'IMPORTANT 2:' below).
# (C) Create a shell script to post-process the results, that will run 
#     after all the jobs have finished running.
#
# IMPORTANT: This script may need to be modified to suit the user's purposes. 
#            See comments below that start with IMPORTANT X, where X is a number 
#            for details.

require "sre.pl";
use Getopt::Std;
use Cwd;
$e_cutoff = 2;
$b_cutoff = 0.0;
$use_evalues   = 1;
$use_bitscores = 0;
$option_m_enabled = 0;
$option_m = "";

getopts('E:B:M:');
if (defined $opt_E) { $e_cutoff = $opt_E; }
if (defined $opt_B) { $b_cutoff = $opt_B; $use_evalues = 0; $use_bitscores = 1; }
$option_m = "";
if (defined $opt_M) { $option_m_enabled = 1; $cm_dir = $opt_M; }

$usage = "Usage: perl rmark_clusterfy.pl\n\t<.rmm file name>\n\t<.rmk file name>\n\t<dir with *.ali *.idx *.test and *.raw files (and *ali.cm if infernal-given-cm.rmm is used)>\n\t<index file with family names; provide path>\n\t<genome root X, X.chrlist, X.fa, X.ebd must be in seq dir>\n\t<output file root>\n";
$options_usage  = "\nOptions:\n\t";
$options_usage .= "-E <x> : use E-values [default], set max E-val to keep as <x> [df: 2]\n\t";
$options_usage .= "-B <x> : use bit scores, set min score to keep as <x>\n\t";
$options_usage .= "-M <x> : pass -M <x> to .rmm module use, to use CM files in dir <x> instead of building them as needed\n\n";

if(@ARGV != 6)
{
    print $usage;
    print $options_usage;
    exit();
}

($rmm, $rmk, $seq_dir, $fam_idx, $genome_root, $out_file_root) = @ARGV;
$seq_dir = getcwd() . "\/" . $seq_dir;
if($option_m_enabled) { 
    $cm_dir = getcwd() . "\/" . $cm_dir; 
    $option_m = "-M $cm_dir";
}

$orig_rmk = $rmk;

# Make a new directory where the benchmark will run ($run_dir)
$run_dir = $out_file_root . "_" . $genome_root . "_out_dir";
if(! (-e "$run_dir/")) { system("mkdir $run_dir" ); } 

# Ensure that files we need are in the seq directory and
# copy them into the dir we're going to run from.
#
# IMPORTANT 2: Non-infernal may not require these files and may require other files, 
#              delete those that are not needed and add new ones here!
#
$chrom_list    = $seq_dir . "/" . $genome_root . ".chrlist";
$genome_file   = $seq_dir . "/" . $genome_root . ".fa";
$embed_file    = $seq_dir . "/" . $genome_root . ".ebd";

if(! (-e ("rmark.pl"))) { die("ERROR, rmark.pl must exist in the current directory."); } 
else { system("cp rmark.pl $run_dir"); } 
if(! (-e ("rmark_process_glbf.pl"))) { die("ERROR, rmark_process_glbf.pl must exist in the current directory."); } 
else { system("cp rmark_process_glbf.pl $run_dir"); } 
if(! (-e ("rmark_times.pl"))) { die("ERROR, rmark_process_glbf.pl must exist in the current directory."); } 
else { system("cp rmark_times.pl $run_dir"); } 
if(! (-e ("infernal.pm"))) { die("ERROR, infernal.pm must exist in the current directory."); } 
else { system("cp infernal.pm $run_dir"); } 
if(! (-e ("infernal2glbf.pl"))) { die("ERROR, infernal2glbf.pl must exist in the current directory."); } 
else { system("cp infernal2glbf.pl $run_dir"); }
if(! (-e ("sre.pl"))) { die("ERROR, sre.pl must exist in the current directory."); } 
else { system("cp sre.pl $run_dir"); } 

if(! ( -e ("$genome_file"))) { die("ERROR, $genome_file must exist in $seq_dir") } 
if(! ( -e ("$embed_file"))) { die("ERROR, $embed_file must exist in $seq_dir.") }
if(! ( -e ("$chrom_list"))) { die("ERROR, $chrom_list must exist in $seq_dir.") }
if(! ( -e ("$fam_idx"))) { die("ERROR, $fam_idx must exist in $seq_dir.") }
if(! (-e ("$rmm"))) { die("ERROR, $rmm doesn't exist."); } 
else { system("cp $rmm $run_dir"); }
if(! (-e ("$rmk"))) { die("ERROR, $rmk doesn't exist."); }
else { system("cp $rmk $run_dir"); }

#IMPORTANT, we've copied the $rmm and $rmk files, now make sure they
# don't include a full path to the file, just the name
$rmm =~ s/.+\///;
$rmk =~ s/.+\///;
$fam_idx_root = $fam_idx;
$fam_idx_root =~ s/.+\///;

#if(! ( -e ("$genome_file"))) { die("ERROR, $genome_file must exist in the current directory.") }
#else { system("cp $genome_file $run_dir"); } 
#if(! ( -e ("$embed_file"))) { die("ERROR, $embed_file must exist in the current directory.") }
#system("cp $embed_file $run_dir");
#if(! ( -e ("$rmm"))) { die("ERROR, $rmm must exist in the current directory.") }
#system("cp $rmm $run_dir");
#if(! ( -e ("$rmk"))) { die("ERROR, $rmk must exist in the current directory.") }
#system("cp $rmk $run_dir");
#if(! ( -e ("$fam_idx"))) { die("ERROR, $fam_idx must exist in the current directory.") }
#if(! ( -e ("$chrom_list"))) { die("ERROR, $chrom_list must exist in the current directory.") }

# Copy any .prior files we might need.
if(-e ("*.pri*")) { system("cp *.pri* $run_dir"); }
# Copy any .null files we might need.
if(-e ("*.null*")) { system("cp *.null* $run_dir"); }
    
# Read in the roots of the chromosomes and determine the size of the genome from chromosome sizes in chromosome list file
file_lines_to_arr($chrom_list, \@chrom_files_arr);
open(IN, $chrom_list) || die("Unable to open chromosome list file $chrom_list\n");
@chrom_files_arr = ();
$total_genome_size = 0;
while($line = <IN>)
{
    chomp $line;
    if($line ne "")
    {
	@els = split(/\s+/, $line);
	if(scalar(@els) != 2) { die("ERROR, line $line of chromosome list file $chrom_list is invalid. Should have exactly two tokens separated by whitespace <chromsome root name> <size of chromosome in nt>\n"); }
	push(@chrom_files_arr, $els[0]);
	$total_genome_size += $els[1];
    }
}
close(IN);
$total_genome_search_space = $total_genome_size * 2.; # we search both strands with cmsearch

# Ensure that if we're using E-values for cmsearch, we've set -Z <dbsize> on the command line, with the
# appropriate dbsize (twice the sum of all the chromosome lengths given in $chrom_list), if not die.
# We do this so that when we split up the chromosomes into one per file, we know that cmsearch will
# report E-values as if the database size was the full pseudogenome, not just the size of one chromosome.
# This makes a clusterfied run, where each search is 1 chromosome, report identical E-values to a non-clusterfied
# search, where each search is against the full genome.
$found_cmsearch = 0;
$found_z_option = 0;
$z_size_agrees  = 0;
if($use_evalues) { 
    open(RMK, $orig_rmk); 
    while($line = <RMK>) { 
	chomp $line;
	if($line =~ /\$cms/) { 
	    $found_cmsearch = 1;
	    if($line =~ /\$cms\s*\=\s*\".+\-Z\s+(\S+).*\"/) { #cmsearch ... -Z <dbsize> ...
		$found_z_option = 1;
		$z_size = $1;
	    }
	    if($line =~ /\$cms\s*\=\s*\".+\-Z\s*\=\s*(\S+).*\"/) { #cmsearch ... -Z\s*=\s*<dbsize> ...
		$found_z_option = 1;
		$z_size = $1;
	    }
	}
	elsif($line !~ m/\#/) { push(@rmk_lines_arr, $line); }
    }
}
if($found_cmsearch) { 
    if(!($found_z_option)) { die("ERROR, using E-values for cmsearch run (found \$cms definition in rmk file\n\t$rmk) but -Z in not enabled with -Z <dbsize>. -Z is required\n\tso when genome is split up into single chromosomes and searched, cmsearch\n\treports E-values as if the entire genome is being searched.\n"); }
    $total_genome_search_space_Mb = ($total_genome_search_space / 1000000.);
    if((abs($z_size - $total_genome_search_space_Mb)) > 0.000001) { 
	die("ERROR, using E-values for cmsearch run (found \$cms definition in rmk file\n\t$rmk) but <dbsize> from rmk file with -Z <dbsize> is invalid.\n\t<dbsize> should be $total_genome_search_space_Mb (2 * the summed chromosome sizes from\n\t$chrom_list\n\t(b/c we're searching both strands)), but read $z_size in rmk file.\n");
    }
}

# Read in the roots of the test families
file_lines_to_arr($fam_idx, \@fam_roots_arr);

# Create the script for the cluster that will submit
# all the jobs. 
push(@exec_lines, "#!/bin/sh");
# For each family...
for($i = 0; $i < scalar(@fam_roots_arr); $i++)
{
    $fam = $fam_roots_arr[$i];
    $num = $i + 1;
    $index_file = "$run_dir/INDEX" . $num;
    # Create a INDEX file for rmark.pl to read.
    open(OUT, ">" . $index_file);
    print OUT ("$fam\n");
    close(OUT);
    # For each chromosome...
    for($j = 0; $j < scalar(@chrom_files_arr); $j++)
    {
	#Determine the chromosome file, and check that it exists in the test_dir.
	$chrom_file = $chrom_files_arr[$j];
	if(! (-e ("$seq_dir" . "\/$chrom_file")))
	{ die("ERROR, $chrom_file must exist in " . $seq_dir. "\n"); }
	$rmark_output = $out_file_root . "_" . $fam_roots_arr[$i] . "_" . $chrom_files_arr[$j];
	$cluster_index_file = $index_file;
	$cluster_index_file =~ s/.+\///;
	$job_name = "rm-$fam-$j";
	# Create the rmark.pl executing line for this family, this chromosome.
	if($use_evalues)
	{
	    $rmark_call = "rmark.pl $option_m -E " . $e_cutoff . " " . $rmm . " " . $rmk . " " . $seq_dir . " " . $cluster_index_file . " " . $chrom_files_arr[$j] . " " . $rmark_output;
	    }
	elsif($use_bitscores)
	{
	    $rmark_call = "rmark.pl $option_m -B " . $b_cutoff . " " . $rmm . " " . $rmk . " " . $seq_dir . " " . $cluster_index_file . " " . $chrom_files_arr[$j] . " " . $rmark_output;
	}
	# Create a command the cluster will make to run rmark.pl
	# IMPORTANT 1: this is a version 6 SGE qsub command - works at Janelia Farm; 
	#              not sure about elsewhere...
	$exec_line = "qsub -q c05.q,c06.q,c07.q,c08.q,c09.q,c10.q,c11.q,c12.q,c13.q,c14.q -N $job_name -o /dev/null -b y -cwd -V -j y \'perl " . $rmark_call . "\'";
	push(@exec_lines, $exec_line);
    }
}	
$command_file_name = $out_file_root . ".com";
print_arr_to_file(\@exec_lines, $command_file_name);
print_out_file_notice($run_dir . "/" . $command_file_name, "Command file with " . (scalar(@exec_lines)-1) . " qsub calls for the cluster.");

# Move the command file into the dir we're going to run from:
system("mv $command_file_name $run_dir");

# Create a shell script to post-process the results to run after
# all the jobs have finished running.
$pp_file = $out_file_root . "_pp.script";
open(PP, ">" . $pp_file);
print PP ("rm merged_" . $out_file_root . "*\n");
$all_glbf_out = $out_file_root . "_all_glbf.concat";
$all_time_out = $out_file_root . "_all_time.concat";
print PP ("cat *.glbf > $all_glbf_out\n");
print PP ("cat *.time > $all_time_out\n");
#11.25.05 - get timing info
print PP ("perl rmark_times.pl *.time > merged_" . $out_file_root . ".time\n");

# Call rmark_process_glbf.pl with defaults: 'hit' resolution mode and 
# ignore cross-hits on both strands.
if($use_evalues)
{
    $rmark_process_option = "E";
}
else
{
    $rmark_process_option = "B";
}
print PP ("perl rmark_process_glbf.pl $rmark_process_option $rmm $rmk $seq_dir $fam_idx_root $genome_root $all_glbf_out merged_" . $out_file_root . "_hit\n");


# Here's some alternative rmark_process_glbf.pl calls :
# the following 2 lines get results in 'fnt' and 'nnt' resolution modes 
#print PP ("perl rmark_process_glbf.pl -R fnt $rmark_process_option $rmm $rmk  $seq_dir $fam_idx_root $genome_root $all_glbf_out merged_" . $out_file_root . "_fnt\n");
#print PP ("perl rmark_process_glbf.pl -R nnt $rmark_process_option $rmm $rmk  $seq_dir $fam_idx_root $genome_root $all_glbf_out merged_" . $out_file_root . "_nnt\n");

#the following rmark_process_glbf.pl calls DON'T IGNORE CROSS HITS!
#print PP ("perl rmark_process_glbf.pl -I none $rmark_process_option $rmm $rmk  $seq_dir $fam_idx_root $genome_root $all_glbf_out merged_" . $out_file_root . "_hit_cross\n");
#print PP ("perl rmark_process_glbf.pl -I none -R fnt $rmark_process_option $rmm $rmk  $seq_dir $fam_idx_root $genome_root $all_glbf_out merged_" . $out_file_root . "_fnt_cross\n");
#print PP ("perl rmark_process_glbf.pl -I none -R fnt nnt $rmark_process_option $rmm $rmk  $seq_dir $fam_idx_root $genome_root $all_glbf_out merged_" . $out_file_root . "_nnt_cross\n");

#the followign rmark_process_glbf.pl calls IGNORE CROSS HITS ONLY ON THE SAME STRAND!
#print PP ("perl rmark_process_glbf.pl -I opp $rmark_process_option $rmm $rmk  $seq_dir $fam_idx_root $genome_root $all_glbf_out merged_" . $out_file_root . "_hit_samecross\n");
#print PP ("perl rmark_process_glbf.pl -I opp -R fnt $rmark_process_option $rmm $rmk  $seq_dir $fam_idx_root $genome_root $all_glbf_out merged_" . $out_file_root . "_fnt_samecross\n");
#print PP ("perl rmark_process_glbf.pl -I opp -R nnt $rmark_process_option $rmm $rmk  $seq_dir $fam_idx_root $genome_root $all_glbf_out merged_" . $out_file_root . "_nnt_samecross\n");

print PP ("cp merged_*fam ../\n");
print PP ("cp merged_*all ../\n");
print PP ("cp merged_*time ../\n");
print PP ("cp merged_*roc ../\n");
close(PP);
print_out_file_notice($run_dir . "/" . $pp_file, "Shell script to merge and process the collective output\n               after all the cluster jobs are finished.");
system("mv $pp_file $run_dir");
system("cp $genome_file $run_dir");
system("cp $fam_idx $run_dir");
system("cp $embed_file $run_dir");

system("chmod +x $run_dir/$pp_file");
# END OF SCRIPT
#################################################################

#################################################################
# Subroutines called in script:
#################################################################
# subroutine : file_lines_to_arr
# from       : M_gen.pm
#
# EPN 03.08.05
#
# purpose : Open a file, and (after chomping) add each line
#           to an array that was passed in.
#
# args (1) $file_in
#          name of file to open and read
#      (2) $arr_ref
#          reference to the array that we'll fill
#################################################################
sub file_lines_to_arr
{
    ($file_in, $arr_ref) = @_;
    open(IN, $file_in) || die;
    while($line = <IN>)
    {
	chomp $line;
	if($line ne "")
	{
	    push(@{$arr_ref}, $line);
	}
    }
}

#################################################################
# subroutine : print_arr_to_file
# from       : M_gen.pm
# 
# EPN 05.25.05
# 
# purpose : Print to a file the elements of a 
#           given array, each on a separate line.
#
# args : (1) $arr_ref 
#            reference to array to print
#        (2) $out_file
#            name of file to print to
################################################################# 
sub print_arr_to_file
{
    ($arr_ref, $out_file)= @_;
    open(OUT, ">" . $out_file);
    for($i = 0; $i < scalar(@{$arr_ref}); $i++)
    { print OUT "$arr_ref->[$i]\n"; }
    close(OUT);
}

#################################################################
# subroutine : print_out_file_notice
# from       : M_gen.pm
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
    $spec_line = "";
    for($i = 0; $i < 75; $i++)
    { $spec_line .= $char; }
    $spec_line .= "\n";
    print("$spec_line");
    print(" Output file notice\n");
    print(" File name   : $file_name\n");
    print(" description : $description\n");
    print("$spec_line");
}


