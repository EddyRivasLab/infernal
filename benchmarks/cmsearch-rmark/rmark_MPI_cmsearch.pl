#!/usr/local/bin/perl
#
# EPN, Mon Dec 11 14:27:41 2006
# rmark_MPI_cmsearch.pl
#
# Prepare a MPI enabled cmsearch RMARK benchmark 
# for running on the cluster.
#
# Usage: perl rmark_MPI_cmsearch.pl
#             <num procs to use>
#             <.rmm file name>
#             <.rmk file name>
#             <seq dir with *.ali *.idx *.test and *.raw files>
#             <index file with fam names; provide path>
#             <genome root X, X.fa, X.ebd must be in seq dir>
#             <dir with CM files, must have <fam>.cm for each <fam>>
#             <output file root>
# Options:
#        -E <x> : use E-values [default], set max E-val to keep as <x> [df: 100]
#        -B <x> : use bit scores, set min score to keep as <x>
#        -O <x> : using old, version 0.x Infernal [default: using 1.x Infernal]
#        -A     : use all nodes of cluster

# Example:  perl rmark_MPI_cmsearch.pl 100 infernal.rmm inf-71.rmk rmark-test/ rmark-test 
#                                      inf-71
#
# The example run above will create the following:
#     - a inf-71_rmark-test_out_dir directory with all the files needed to run
#       rmark.pl copied to it (except seq files which stay in <seq dir>)
#
# General Strategy 
# (A) Create a script to submit X jobs (one for each family) to the 
#     cluster to run cmsearch in MPI mode. 
# (B) Copy all required files into a new directory where the cmsearch jobs 
#     will be run (SEE 'IMPORTANT 2:' below).
# (C) Create a shell script to post-process the results, that will be run 
#     locally after all the jobs have finished running.
#
# IMPORTANT: This script may need to be modified to suit the user's purposes. 
#            See comments below that start with IMPORTANT X, where X is a number 
#            for details.

require "sre.pl";
use Getopt::Std;
use Cwd;
$e_cutoff = 100;
$b_cutoff = 0.0;
$use_evalues   = 1;
$use_bitscores = 0;
$pre_version1 = 0;
$use_all_nodes = 0;

getopts('E:B:O:A');
if (defined $opt_E) { $e_cutoff = $opt_E; }
if (defined $opt_B) { $b_cutoff = $opt_B; $use_evalues = 0; $use_bitscores = 1; }
if (defined $opt_O) { $pre_version1 = 1; }
if (defined $opt_A) { $use_all_nodes = 1; }

$usage = "Usage: perl rmark_MPI_cmsearch.pl\n\t<num processors to use>\n\t<.rmm file name>\n\t<.rmk file name>\n\t<dir with *.ali *.idx *.test and *.raw files>\n\t<index file with family names; provide path>\n\t<genome root X, X.fa, X.ebd must be in seq dir>\n\t<directory with CM files>\n\t<output file root>\n";
$options_usage  = "\nOptions:\n\t";
$options_usage .= "-E <x> : use E-values [default], set max E-val to keep as <x> [df: 2]\n\t";
$options_usage .= "-B <x> : use bit scores, set min score to keep as <x>\n\t";
$options_usage .= "-O <x> : using old, version 0.x Infernal [default: using 1.x Infernal]\n\n";
$options_usage .= "-A     : use all nodes of cluster\n\n";

if(@ARGV != 8)
{
    print $usage;
    print $options_usage;
    exit();
}

($nprocs, $rmm, $rmk, $seq_dir, $fam_idx, $genome_root, $cm_dir, $out_file_root) = @ARGV;
$seq_dir = getcwd() . "\/" . $seq_dir;
$cm_dir  = getcwd() . "\/" . $cm_dir;
$orig_rmk = getcwd() . "\/" . $rmk;

# Make a new directory where the benchmark will run ($run_dir)
$run_dir = $out_file_root . "_" . $genome_root . "_out_dir";
if(! (-e "$run_dir/")) { system("mkdir $run_dir" ); } 

# Ensure that files we need are in the seq directory and
# copy them into the dir we're going to run from.
#
# IMPORTANT 2: Non-infernal may not require these files and may require other files, 
#              delete those that are not needed and add new ones here!
#
$genome_file   = $seq_dir . "/" . $genome_root . ".fa";
$embed_file    = $seq_dir . "/" . $genome_root . ".ebd";

if(! (-e ("rmark_process_glbf.pl"))) { die("ERROR, rmark_process_glbf.pl must exist in the current directory."); } 
else { system("cp rmark_process_glbf.pl $run_dir"); } 
if(! (-e ("rmark_times.pl"))) { die("ERROR, rmark_times.pl must exist in the current directory."); } 
else { system("cp rmark_times.pl $run_dir"); } 
if(! (-e ("infernal2time.pl"))) { die("ERROR, infernal2time.pl must exist in the current directory."); } 
else { system("cp infernal2time.pl $run_dir"); }
if(! (-e ("infernal.pm"))) { die("ERROR, infernal.pm must exist in the current directory."); } 
else { system("cp infernal.pm $run_dir"); } 
if(! (-e ("infernal2glbf.pl"))) { die("ERROR, infernal2glbf.pl must exist in the current directory."); } 
else { system("cp infernal2glbf.pl $run_dir"); }
if(! (-e ("sre.pl"))) { die("ERROR, sre.pl must exist in the current directory."); } 
else { system("cp sre.pl $run_dir"); } 

if(! ( -e ("$genome_file"))) { die("ERROR, $genome_file must exist in $seq_dir") } 
if(! ( -e ("$embed_file"))) { die("ERROR, $embed_file must exist in $seq_dir.") }
if(! ( -e ("$fam_idx"))) { die("ERROR, $fam_idx must exist in $seq_dir.") }
if(! (-e ("$rmm"))) { die("ERROR, $rmm doesn't exist."); } 
else { system("cp $rmm $run_dir"); }
if(! (-e ("$rmk"))) { die("ERROR, $rmk doesn't exist."); }
else { system("cp $rmk $run_dir"); }

require("$orig_rmk");
if($cms eq "")
{
    die("ERROR, in MPI mode the RMARK config file $orig_rmk\nmust define variable \$cm\n");
}

#IMPORTANT, we've copied the $rmm and $rmk files, now make sure they
# don't include a full path to the file, just the name
$rmm =~ s/.+\///;
$rmk =~ s/.+\///;
$fam_idx_root = $fam_idx;
$fam_idx_root =~ s/.+\///;

# Copy any .prior files we might need.
system("cp *.pri* $run_dir");
# Copy any .null files we might need.
system("cp *.null* $run_dir"); 

# Read in the roots of the test families
file_lines_to_arr($fam_idx, \@fam_roots_arr);

# For each family, build a CM using specs in 
# rmk file. THIS IS NO LONGER DONE, WE REQUIRE
# THE CM ALREADY BUILT. CODE LEFT HERE FOR REFERENCE
#for($i = 0; $i < scalar(@fam_roots_arr); $i++)
#{
#    $fam = $fam_roots_arr[$i];
#   $num = $i + 1;
#    $cm_file = $fam . ".cm";
#    $ali_file = $seq_dir . "/" . $fam . ".ali";
#    system("$cmb $fam.cm $ali_file > /dev/null");
#    system("mv $fam.cm $run_dir");
#    printf("mv $fam.cm $run_dir\n");
#}

# Create the script for the cluster that will submit
# the cmsearch jobs, one job for each family.
push(@exec_lines, "#!/bin/sh");
# For each family...
for($i = 0; $i < scalar(@fam_roots_arr); $i++)
{
    $fam = $fam_roots_arr[$i];
    $cm  = $cm_dir . "/" . $fam . ".cm";
    if(! (-e "$cm")) { die("ERROR cm $cm does not exist.\n"); } 
    $job_name = "rm-$fam";
    $cmsearch_name = $run_dir . "\/" . "rm-$fam.cmsearch";
    $out_name = $run_dir . "\/" . "rm-$fam.out";

    if($pre_version1) 
    { 
	#$cmsearch_call = "mpirun -l C $cms --noalign $cm $genome_file"; 
	$cmsearch_call = "mpirun -mca btl self,tcp --prefix /usr/local/openmpi -np $nprocs $cms --noalign $cm $genome_file"; 
    }
    else 
    { 
	#$cmsearch_call = "mpirun -l C $cms --mpi --noalign $cm $genome_file"; 
	$cmsearch_call = "mpirun -mca btl self,tcp --prefix /usr/local/openmpi -np $nprocs $cms --mpi --noalign $cm $genome_file"; 
    }

    if($use_all_nodes) 
    {
	#$exec_line = "qsub -N $job_name -o $out_name -b y -cwd -V -j y -pe lam-mpi-tight $nprocs \'" . $cmsearch_call . " > $cmsearch_name\'";
	$exec_line = "qsub -N $job_name -o $out_name -b y -cwd -V -j y -pe openmpi $nprocs \'" . $cmsearch_call . " > $cmsearch_name\'";
    }
    else # only use c05-c14 
    {
	#$exec_line = "qsub -q c05.q,c06.q,c07.q,c08.q,c09.q,c10.q,c11.q,c12.q,c13.q,c14.q -N $job_name -o $out_name -b y -cwd -V -j y -pe lam-mpi-tight $nprocs \'" . $cmsearch_call . " > $cmsearch_name\'";
	$exec_line = "qsub -q c05.q,c06.q,c07.q,c08.q,c09.q,c10.q,c11.q,c12.q,c13.q,c14.q -N $job_name -o $out_name -b y -cwd -V -j y -pe openmpi $nprocs \'" . $cmsearch_call . " > $cmsearch_name\'";
    }
    push(@exec_lines, $exec_line);
}


$command_file_name = getcwd() . "\/" . $out_file_root . ".com";
print_arr_to_file(\@exec_lines, $command_file_name);
print_out_file_notice($command_file_name, "Command file with " . (scalar(@exec_lines)-1) . " qsub calls for the cluster.");

# 01.22.07 Keep the command file where it is, dammit.
# Move the command file into the dir we're going to run from:
#system("mv $command_file_name $run_dir");
system("cp $command_file_name $run_dir");
    
# Create a shell script to post-process the results to run after
# all the jobs have finished running.
$pp_file = $out_file_root . "_pp.script";
open(PP, ">" . $pp_file);
# First we need to get the glbf files:
$time_name = $run_dir . "\/" . "rm_time.concat";
print PP ("grep \"time\" " . $run_dir . "\/*.cmsearch > $time_name\n");
$all_time_out = $out_file_root . ".time";
print PP ("perl infernal2time.pl $time_name > $run_dir/$all_time_out\n");
print PP ("cp $run_dir/$all_time_out ./\n");

for($i = 0; $i < scalar(@fam_roots_arr); $i++)
{
    $fam = $fam_roots_arr[$i];
    $cmsearch_name = $run_dir . "\/" . "rm-$fam.cmsearch";
    $glbf_name = $run_dir . "\/" . "rm-$fam.glbf";

    if($use_evalues)
    {
	print PP ("echo '>$fam' > $glbf_name\n");
	print PP ("perl infernal2glbf.pl -E $e_cutoff $cmsearch_name >> $glbf_name\n");
    }
    else
    {
	print PP ("echo '>$fam' > $glbf_name\n");
	print PP ("perl infernal2glbf.pl -B $b_cutoff $cmsearch_name >> $glbf_name\n");
    }
}

# Now we build the PP (post-processing) script the same as in
# rmark_clusterfy.pl.
print PP ("rm merged_" . $out_file_root . "_hit*\n");
$all_glbf_out = $out_file_root . "_all_glbf.concat";
#$all_time_out = $out_file_root . "_all_time.concat";
print PP ("cat "  . $run_dir . "\/" . "*.glbf > $all_glbf_out\n");
#print PP ("cat "  . $run_dir . "\/" . "*.time > $all_time_out\n");
#11.25.05 - get timing info
#print PP ("perl rmark_times.pl " . $run_dir . "\/" . "*.time > merged_" . $out_file_root . ".time\n");

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

#print PP ("cp merged_*fam ../\n");
#print PP ("cp merged_*all ../\n");
#print PP ("cp merged_*time ../\n");
#print PP ("cp merged_*roc ../\n");
close(PP);
print_out_file_notice($pp_file, "Shell script to merge and process the collective output\n               after all the cluster jobs are finished.");
# 01.22.07 Keep the pp file where it is, dammit.
#system("mv $pp_file $run_dir");
system("cp $pp_file $run_dir");
system("cp $genome_file $run_dir");
system("cp $fam_idx $run_dir");
system("cp $embed_file $run_dir");

system("chmod +x $pp_file");
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


