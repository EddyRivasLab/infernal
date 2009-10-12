#! /usr/local/bin/perl
#
# Eric Nawrocki 09.26.05
# rmark.pl 
#
# Runs a RMARK benchmark, using the search software defined in the
# .rmm module and options defined in the .rmk file.
#
# Usage:    perl rmark.pl 
#                <.rmm module> 
#                <.rmk config file>
#                <seq directory>
#                <index file with fam names; provide path>
#                <genome file; must be seq dir>
#                <output root, for naming output files>
#
# Options:
#        -E <x> : use E-values [default], set max E-val to keep as <x> [default: 2]
#        -B <x> : use bit scores, set min score to keep as <x>
#        -M <x> : pass -M <x> to .rmm module use, to use CM files in dir <x> instead of building them as needed
#        -S <x> : (with rfamscan.rmm) save temporary files with root name <x>
#
# Example:  perl rmark.pl infernal.rmm inf-72.rmk rmark-test/ rmark-test.idx
#                         rmark-test.fa rmark-test_out
#
# The example run above will result in the following files:
# rmark_test_out.glbf: the glbf output derived from the blast/infernal output.
# rmark_test_out.time: the running time of the program.
#
# The *.glbf file can be transformed into results for several different
# scoring schemes by rmark_process_glbf.pl
#
# For each family name fam listed in X.idx, there must be these files:
#      <seq directory>/fam.idx  : list of true positive test seqs
#      <seq directory>/fam.ali  : alignment of training seqs
#      <seq directory>/fam.test : FASTA file of test seqs
#      <seq directory>/fam.raw  : FASTA file of unaligned (raw) training seqs
#
use Getopt::Std;
$e_cutoff = 1;
$b_cutoff = 0.0;
$use_evalues   = 1;
$use_bitscores = 0;

getopts('E:B:M:S:');
if (defined $opt_E) { $e_cutoff = $opt_E; }
if (defined $opt_B) { $b_cutoff = $opt_B; $use_evalues = 0; $use_bitscores = 1; }
$option_m = "";
$option_s = "";
if (defined $opt_M) { $option_m = "-M $opt_M"; }
if (defined $opt_S) { $option_s = "-S $opt_S"; }

$usage = "Usage: perl rmark.pl\n\t<.rmm rmark module>\n\t<.rmk rmark config file>\n\t<seq directory with *.ali, *.test, *.idx, *.raw files>\n\t<index file with family names; provide path>\n\t<genome file; must be in seq dir>\n\t<output root, for naming output files>\n";
$options_usage  = "\nOptions:\n\t";
$options_usage .= "-E <x> : use E-values [default], set max E-val to keep as <x> [default: 1]\n\t";
$options_usage .= "-B <x> : use bit scores, set min score to keep as <x>\n\t";
$options_usage .= "-M <x> : pass -M <x> to .rmm module use, to use CM files in dir <x> instead of building them as needed\n\t";
$options_usage .= "-S <x> : (with rfamscan.rmm) save temporary files with root name <x>\n\n";

if(@ARGV != 6)
{
    print $usage;
    print $options_usage;
    exit();
}

$rmm = shift;
$rmk = shift;
$dir = shift;
$idx = shift;
$genome_file = shift;
$out_root = shift;
$genome_file = $dir . "/" . $genome_file;
$total_runtime = 0;

open(GLBF, ">" . $out_root . ".glbf");
open(TIME, ">" . $out_root . ".time");

($full_bsec, $full_bmin, $full_bhour, $full_bdate, $full_bmonth, $full_byear, $full_bweekday, $full_byearday, $full_bisdst) = localtime;
open (INDEX,$idx) || die;
while (<INDEX>) {
    if (/^(\S+)/) {
	$fam = $1;
	printf GLBF (">$fam\n");

	$cur_option_m = "";
	if($option_m ne "") { 
	    if($option_m =~ m/\/$/) { $cur_option_m = $option_m . "$fam.cm";  }
	    else                    { $cur_option_m = $option_m . "/$fam.cm"; }
	}
	# Run the search module
	($bsec, $bmin, $bhour, $bdate, $bmonth, $byear, $bweekday, $byearday, $bisdst) = localtime;
	if($use_evalues)
	{
	    printf("perl $rmm $option_s $cur_option_m -E $e_cutoff $rmk $dir/$fam.idx $dir/$fam.ali $genome_file\n");
	    $glbfoutput = `perl $rmm $cur_option_m -E $e_cutoff $rmk $dir/$fam.idx $dir/$fam.ali $genome_file`;
	}
	elsif($use_bitscores)
	{
	    $glbfoutput = `perl $rmm $cur_option_m -B $b_cutoff $rmk $dir/$fam.idx $dir/$fam.ali $genome_file`;
	}    
	($esec, $emin, $ehour, $edate, $emonth, $eyear, $eweekday, $eyearday, $eisdst) = localtime;
	# We calculate run time as the elapsed time with
	# resolution only at the seconds level. This isn't robust if the month changes during
        # execution.
	$ddate = $edate - $bdate;
	$dhour = $ehour - $bhour;
	$dmin =  $emin - $bmin;
	$dsec =  $esec - $bsec;
	$sec_runtime = $ddate * 24 * 60 * 60;
	$sec_runtime += $dhour * 60 * 60;
	$sec_runtime += $dmin * 60;
	$sec_runtime += $dsec;
	
	#print("glbfoutput:\n$glbfoutput\n");

	# GLBF output will have all hits to each sequence
        # GLBF format is just like GLF format but with bounds of hits
	#      and with orientation of hits (0 for forward strand, 1 for reverse)
        # <seq name> <score> <start posn> <end posn> <orientation>
	
	@lines = split(/^/, $glbfoutput);
	foreach $line (@lines) 
	{
	    print GLBF $line;
	}
	printf TIME ("$fam:search_runtime(secs_elapsed): $sec_runtime\n");
	$total_runtime += $sec_runtime;
    }
}
($full_esec, $full_emin, $full_ehour, $full_edate, $full_emonth, $full_eyear, $full_eweekday, $full_eyearday, $full_eisdst) = localtime;

printf TIME ("total:search_runtime(secs_elapsed): $total_runtime\n");

close(GLBF);
close(TIME);

